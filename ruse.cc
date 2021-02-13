/* ruse.cc -- modify a file without changing its CRC

  Copyright (C) 2018, 2021 Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Mark Adler
  madler@alumni.caltech.edu

 */

  /*
   Change a few bits in a file such that the CRC of the file doesn't change.

   Usage:

       ruse deg poly ref path

   where deg is the degree of the CRC polynomial in decimal (i.e. the length of
   the CRC in bits); poly is the CRC polynomial in hexadecimal, not including
   the x^deg term; ref is 0 or 1, where 1 indicates that the CRC is reflected;
   and path is the path of the file to be modified. For example:

       ruse 32 4c11db7 1 file

   will modify file so that the standard CRC-32 of file is the same before and
   after. ruse will find the codeword it knows that has the fewest changes that
   it can apply, often as few as three bit flips for files of tens of Kbytes or
   more. It will then apply the codeword starting at a randomly chosen bit
   location in the file. ruse loads a few thousand codewords from the file
   codewords.txt, which were scraped from Philip Koopman's web site on CRC
   codewords. (See the associated getcodes script.) If the provided polynomial
   has no codewords from that site, or if the codewords.txt file is not found,
   then the polynomial itself is used as the codeword.

   If no arguments are given, then the number of CRCs and codewords loaded from
   codewords.txt is shown.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <regex>
#include <random>
#include <limits>
#include <cctype>
#include <cstdint>
#include <climits>
using namespace std;
#define local static

// Any carets (^) appearing in the comments refer to the power operation, as
// opposed to exclusive-or.

// --- bit operations ---

// Return the number of one bits in x. The non-builtin version is efficient for
// sparse ones.
local int weight(uint64_t x) {
#if defined(__has_builtin) && __has_builtin(__builtin_popcountll) && \
    defined(ULLONG_MAX) && (ULLONG_MAX >> 32) == 4294967295UL
    return __builtin_popcountll(x);
#else
    int n = 0;
    while (x) {
        n++;
        x &= x - 1;
    }
    return n;
#endif
}

// Return the bit position of the highest one bit in n (the floor of the log
// base 2 of n), or -1 if n is zero.
local int ilog2(uint64_t n) {
#if defined(__has_builtin) && __has_builtin(__builtin_clzll) && \
    defined(ULLONG_MAX) && (ULLONG_MAX >> 32) == 4294967295UL
    return n ? 63 - __builtin_clzll(n) : -1;
#else
    uint64_t m;
    int k = -1;
    if ((m = n >> 32) != 0) k += 32, n = m;
    if ((m = n >> 16) != 0) k += 16, n = m;
    if ((m = n >> 8) != 0) k += 8, n = m;
    if ((m = n >> 4) != 0) k += 4, n = m;
    if ((m = n >> 2) != 0) k += 2, n = m;
    if ((m = n >> 1) != 0) k++, n = m;
    return k + (int)n;
#endif
}

// Reverse the low n bits of x, setting the remaining bits to zero. The result
// is undefined if n is not in 1..64.
local uint64_t reverse(uint64_t x, int n) {
#if defined(__has_builtin) && __has_builtin(__builtin_bitreverse64)
    return __builtin_bitreverse64(x) >> (64 - n);
#else
    x = (((x & 0xaaaaaaaaaaaaaaaa) >> 1) |  ((x & 0x5555555555555555) << 1));
    x = (((x & 0xcccccccccccccccc) >> 2) |  ((x & 0x3333333333333333) << 2));
    x = (((x & 0xf0f0f0f0f0f0f0f0) >> 4) |  ((x & 0x0f0f0f0f0f0f0f0f) << 4));
    x = (((x & 0xff00ff00ff00ff00) >> 8) |  ((x & 0x00ff00ff00ff00ff) << 8));
    x = (((x & 0xffff0000ffff0000) >> 16) | ((x & 0x0000ffff0000ffff) << 16));
    return ((x >> 32) | (x << 32)) >> (64 - n);
#endif
}

// --- CRC codeword operations ---

// Type for CRC representation and processing. word_t must be an unsigned type
// with a number of bits greater than or equal to the longest CRC to be
// processed. If this type has more than 64 bits, then the bit operations above
// will need to be modified to handle more than 64 bits.
typedef uint64_t word_t;

// This class provides the operations needed to rapidly compute the CRC of a
// message sparse in ones. The class is initialized with the CRC polynomial
// using the Koopman notation (see comments below). The class provides a CRC
// register, crc, and two methods that operate on crc. They are zeros(n), which
// applies n zeros to the CRC, and one(), which applies a single one to the
// CRC. zeros() is written to be very efficient, since the messages of
// interest are sparse, with long series of zeros between occasional ones.
// Sometimes billions of zeros. zeros(n) takes O(log n) time for large n.
struct crc_sparse {
    int deg;                    // degree of the polynomial (bits in CRC)
    word_t poly;                // polynomial, not including x^deg
    word_t rpoly;               // poly reflected
    word_t crc;                 // CRC register value

    // Initialize the CRC description from a polynomial in the Koopman
    // convention. In that convention the high bit of the polynomial (the x^n
    // term for an n-bit CRC) is retained, but the low bit of the polynomial,
    // which must be a 1 (the x^0 term) is down-shifted off the end. This
    // allows a single value that provides the degree of the polynomial, by
    // virtue of the location of the most significant one bit, and that also
    // allows an n-bit CRC description to fit in an n-bit word. Here that
    // convention is converted into what is needed in a generalized CRC
    // implementation, which is the degree of the polynomial as an integer,
    // deg, and the low bits of the polynomial, i.e. including the x^0 term,
    // but not including the x^n term as poly. In addition the bit reverse of
    // poly is generated, rpoly, for reflected implementations of the same CRC.
    // It is very common for software implementations of a CRC to be reflected,
    // since it can save a few operations in the computation.
    //
    // The CRC register value crc is initialized to zero when this class is
    // instantiated.
    //
    // This CRC description is not complete, in that it does not include the
    // initialization of the register, the exclusive-or value applied to the
    // register at the end of a CRC calculation, nor whether the CRC should be
    // reflected at the end of the calculation. Those can be done outside of
    // this class if needed.
    crc_sparse(word_t k) {
        deg = ilog2(k);
        poly = ((k ^ ((word_t)1 << deg)) << 1) | 1;
        deg++;
        rpoly = reverse(poly, deg);
        crc = 0;
    }

  private:
    // Table of matrix operators that apply power-of-two zeros to the CRC.
    // zero_op[n] multiplied by the CRC applies 2^n zeros. This table is built
    // up on demand by zeros().
    vector<vector<word_t>> zero_op;

    // Multiply matrix mat times vector vec over GF(2), returning a vector.
    static word_t mult(vector<word_t> mat, word_t vec) {
        word_t sum = 0;
        unsigned i = 0;
        while (vec) {
            if (vec & 1)
                sum ^= mat[i];
            vec >>= 1;
            i++;
        }
        return sum;
    }

    // Return the square of the matrix mat over GF(2).
    static vector<word_t> square(vector<word_t> mat) {
        vector<word_t> sq;
        sq.reserve(mat.size());
        for (auto vec : mat)
            sq.push_back(mult(mat, vec));
        return sq;
    }

  public:
#   define ZORDERN 12   // number of low bits of n for which to use the O(n)
                        // approach for zeros() -- this parameter is determined
                        // empirically

    // Apply n zeros using the reflected polynomial. This takes O(log n) time
    // for large n, using matrix operators. The operators are generated on
    // demand and retained for later use in zero_op.
    void zeros(uintmax_t n) {
        // The O(n) approach is faster for small n.
        auto k = n & ((1 << ZORDERN) - 1);
        n >>= ZORDERN;
        while (k) {
            crc = crc & 1 ? (crc >> 1) ^ rpoly : crc >> 1;
            k--;
        }

        // Use matrix multiplication for large n. Add entries to zero_op as
        // needed, reusing them in subsequent calls.
        unsigned i = ZORDERN;
        while (n) {
            if (zero_op.size() == 0) {
                // generate operator to apply one zero (2^0 zeros)
                vector<word_t> op;
                op.push_back(rpoly);
                word_t row = 1;
                for (int j = 1; j < deg; j++) {
                    op.push_back(row);
                    row <<= 1;
                }
                zero_op.push_back(op);
            }
            while (zero_op.size() <= i)
                // generate 2^i zeros operator from 2^(i-1) zeros operator
                zero_op.push_back(square(zero_op.back()));
            if (n & 1)
                crc = mult(zero_op[i], crc);
            n >>= 1;
            i++;
        }
    }

    // Apply one one to crc using the reflected polynomial.
    void one() {
        crc = crc & 1 ? crc >> 1 : (crc >> 1) ^ rpoly;
    }
};

// A CRC codeword, in reflected order for application to a message.
struct codeword {
    vector<uintmax_t> ones;     // offsets of the 1 bits in reflected codeword

    // reflect the provided codeword and save that in ones
    codeword(vector<uintmax_t> init) {
        auto n = init.size();
        ones.resize(n);
        auto last = init.back();
        for (auto pos : init)
            ones[--n] = last - pos;
    }
};

// Extend crc_sparse to have a set of codewords associated with the CRC.
struct codewords : public crc_sparse {
    map <uintmax_t, codeword> codes;    // codewords for this CRC (the key is
                                        // the length of the codeword)

    // Add a codeword with ones in the positions in ones[], but only if it
    // really is a codeword. Return true if it is not a codeword.
    int add(vector<uintmax_t> ones) {
        // create the codeword, reflecting what was provided
        auto code = codeword(ones);

        // it is a codeword if and only if it takes a CRC of 0 back to 0
        crc = 0;
        uintmax_t off = 0;
        for (auto pos : code.ones) {
            zeros(pos - off);
            one();
            off = pos + 1;
        }
        if (crc != 0)
            return 1;

        // add the verified codeword to the map, using the length of the
        // codeword in bits as the key -- this choice of key allows for a fast
        // search for the longest codeword less than or equal to a given length
        codes.insert(make_pair(ones.back() + 1, code));
        return 0;
    }

    // Initialize this CRC with a Koopman convention polynomial. Add the first
    // codeword, which is the CRC polynomial itself.
    codewords(word_t k) : crc_sparse(k) {
        // add the polynomial as a codeword
        vector<uintmax_t> ones = {0};           // x^0 term is implied
        while (k) {
            word_t next = k;
            k &= k - 1;
            next ^= k;
            ones.push_back(ilog2(next) + 1);
        }
        if (add(ones))
            throw logic_error("polynomial is not a codeword?!");
    }
};

// Multiply a by b, putting the product in *c. If the multiplication overflows,
// then return true. a, b, and c should not have side effects as they may be
// evaluated more than once.
#if __has_builtin(__builtin_mul_overflow)
#  define MUL(a,b,c) __builtin_mul_overflow(a,b,c)
#else
#  define MUL(a,b,c) (*(c) = (a) * (b), (b) && *(c) / (b) != (a))
#endif

// Convert the string str to a generic integer of type T. Any leading white
// space is skipped, and the number may optionally be preceded by a "+" or "-"
// to indicate the sign. A "-" is allowed for unsigned types, which will result
// in the unary minus operation to be applied to the returned unsigned value.
//
// base must be in the range 0..36, where 1..36 is the base for the digits, or
// 0 indicates that the base is 8, 10, or 16 as determined by the prefix. For
// base == 0, if the first two characters after any leading white space and
// sign is "0x" or "0X", then the base is taken to be hexadecimal (16). If
// there is a leading "0" that is not followed by an "x" or "X", then the base
// is taken to be octal (8). If there is no leading "0", then the base is taken
// to be decimal (10). A leading "0x" or "0X" is also permitted, but not
// required, when the base is given as 16. Note that if the base is 1, then
// only the value 0 can be returned. If base is not provided, base is taken to
// be 10.
//
// If T is not an integer type or base is not in 0..36, or if there are no
// valid digits in the string (where a "0x" or "0X" does not count as a digit),
// then a std::invalid_argument exception is thrown. A single "0" indicating
// octal for base == 0 is considered a digit.
//
// If the number is out of range for the integer type, then a std::out_of_range
// exception is thrown. For an unsigned type, the full range is permitted
// regardless of the sign. For a signed type, only the range of the signed
// type, taking into account the sign, is permitted.
//
// The resulting converted integer of type T is returned. On return, if pos is
// not zero, *pos is set to the index of the first character in the string not
// converted, or the length of the string if all of it was used. If pos is not
// provided, it is taken to be zero.
template<typename T> T stoint(string str, size_t *pos = 0, int base = 10) {
    if (!numeric_limits<T>::is_integer)
        throw invalid_argument("stoint: not an integer type");
    if (base < 0 || base > 36)
        throw invalid_argument("stoint: base out of range");

    // skip any leading white space
    size_t len = str.size();
    size_t i = 0;
    while (i < len && isspace(str[i]))
        i++;

    // process possible sign
    int negate = 0;
    if (i < len) {
        int ch = str[i];
        if (ch == '+')
            i++;
        else if (ch == '-') {
            negate = 1;
            i++;
        }
    }

    // process leading 0 or 0x for octal or hex and resolve base
    int none = 1;                       // true if no digits processed
    if (i < len) {
        if (str[i] == '0') {
            i++;
            if ((base == 0 || base == 16) &&
                i < len && (str[i] == 'x' || str[i] == 'X')) {
                base = 16;
                i++;
            }
            else {
                if (base == 0)
                    base = 8;
                none = 0;               // allow just "0"
            }
        }
    }
    if (base == 0)
        base = 10;

    // process digits
    uintmax_t max = numeric_limits<T>::is_signed && negate ?
                        -(intmax_t)numeric_limits<T>::min() :
                        numeric_limits<T>::max();   // max value after sign
    uintmax_t val = 0;                              // accumlated value
    int u9 = '9' + (base > 10 ? 10 : base) - 10,    // digit upper limits
        uZ = 'Z' + base - 36,
        uz = 'z' + base - 36;
    while (i < len) {
        int dig = str[i];
        dig = dig >= '0' && dig <= u9 ? dig - '0' :
              dig >= 'A' && dig <= uZ ? dig - 'A' + 10 :
              dig >= 'a' && dig <= uz ? dig - 'a' + 10 :
              -1;
        if (dig == -1)
            break;
        uintmax_t mul;
        if (MUL(val, base, &mul) || (val = mul + dig) < mul || val > max)
            throw out_of_range("stoint: overflow");
        none = 0;
        i++;
    }
    if (none)
        throw invalid_argument("stoint: no conversion performed");

    // return position left at in the string and the value as an integer type T
    if (pos)
        *pos = i;
    return negate ? -(T)val : (T)val;
}

// Within this class a catalog of CRCs and their associated codewords is built.
// Methods are provided to add a CRC and its codeword using a specification in
// a provided text string, and to load all of the codewords from a specified
// file. In both cases all provided codewords are verified to be codewords for
// the corresponding CRC polynomials. This class also provides a method to
// modify a file by exclusive-oring it with a codeword of the specified
// polynomial. In so doing, the file will have the same CRC of that polynomial,
// before and after the modification.
struct catalog {
    map<word_t, codewords> cat;

    // Add a codeword from str, which is a line from a Koopman *_len.txt file.
    // Here is a link to one of those files:
    //
    //     https://users.ece.cmu.edu/~koopman/crc/c32/0x82608edb_len.txt
    //
    // Example content from another one, .../crc/c05/0x1e_len.txt:
    // # 0x1e  HD=5  len=1  Example: Len=2 0x1e {0,1} (0x11) (Bits=4)
    // # 0x1e  HD=6  NONE  Example: Len=1 0x1e {0} (0x1e) (Bits=5)
    //
    // In some newer ones, the repeated CRC has been omitted:
    // # 0xa2572962  HD=16  NONE  Example: Len=1 {0} (0xa2572962) (Bits=15)
    //
    // Return true if the str is not in the expected format, or if the data in
    // str is not valid, including if the claimed codeword is not actually a
    // codeword for the provided CRC polynomial.
    int add(string str) {
        // extract the tokens from the string
        istringstream line(str);
        string tok[11];
        for (int i = 0; i < 11; i++)
            line >> tok[i];

        // verify the expected format of the codeword description
        smatch match;               // (not used)
        if (tok[0] != "#" ||
            !regex_match(tok[1], match, regex("0x[0-9a-fA-F]+")) ||
            !regex_match(tok[2], match, regex("HD=\\d+")) ||
            !regex_match(tok[3], match, regex("len=\\d+|NONE")) ||
            tok[4] != "Example:" ||
            !regex_match(tok[5], match, regex("Len=\\d+")))
            return 1;
        if (tok[6] == tok[1]) {
            // skip repeated polynomial if present
            tok[6] = tok[7];
            tok[7] = tok[8];
            tok[8] = tok[9];
            tok[9] = tok[10];
            tok[10] = "";
        }
        if (!regex_match(tok[6], match, regex("\\{\\d+(,\\d+)*\\}")) ||
            !regex_match(tok[7], match, regex("\\(0x[0-9a-fA-F]+\\)")) ||
            !regex_match(tok[8], match, regex("\\(Bits=\\d+\\)")) ||
            tok[9] != "" || tok[10] != "")
            return 1;

        try {
            // extract the relevant contents
            auto k = stoint<word_t>(tok[1], 0, 0);      // Koopman polynomial
            if (k == 0)
                return 1;
            auto len = stoint<uintmax_t>(tok[5].substr(4)); // dataword length
            auto res = stoint<word_t>(tok[7].substr(1), 0, 0);  // residual CRC
            auto bits = stoint<int>(tok[8].substr(6));      // number of ones
            vector<uintmax_t> ones;
            ones.reserve(bits);
            size_t pos = 0;
            for (;;) {
                tok[6] = tok[6].substr(pos + 1);
                if (tok[6].size() == 0)
                    break;
                ones.push_back(stoint<uintmax_t>(tok[6], &pos));    // offset
            }

            // verify the weight
            if (bits < 2 || (int)ones.size() + weight(res) != bits)
                return 1;

            // fill out the rest of the one positions from res
            while (res) {
                word_t next = res;
                res &= res - 1;
                next ^= res;
                ones.push_back(len + ilog2(next));
            }
            if (ones[0] != 0)
                return 1;

            // add the CRC to the catalog, if not there already
            if (cat.find(k) == cat.end())
                cat.insert(make_pair(k, codewords(k)));
            auto crc = &cat.at(k);

            // add the codeword -- return true if it doesn't verify
            return crc->add(ones);
        }
        catch (out_of_range const& e) {
            return 1;
        }
    }

    // Load all of the codewords from path into the catalog. Return the number
    // of invalid lines encountered. If path cannot be opened to read, -1 is
    // returned.
    int load(string path) {
        ifstream in(path);
        if (!in.is_open())
            return -1;
        int bad = 0;
        string str;
        while (getline(in, str))
            if (add(str) && str.substr(0, 2) == "# ")
                bad++;
        return bad;
    }

    // Return the number of CRCs and the total number of codewords in the
    // catalog, in that order.
    pair <int, int> counts() {
        int crcs = 0, codes = 0;
        for (auto const& ent : cat) {
            crcs++;
            codes += ent.second.codes.size();
        }
        return make_pair(crcs, codes);
    }

    // Exclusive-or the longest codeword that fits for the given CRC on the
    // provided file, at a randomly chosen bit offset. The CRC is provided
    // using the Koopman convention, ref is true if the CRC is reflected, and
    // path is the name of the file to modify. If the CRC is not in the
    // catalog, it is added and the CRC polynomial serves as the codeword. The
    // modified file will have the same CRC is the file before modification.
    // The number of bits flipped is returned, or zero if the file is shorter
    // than the shortest codeword, or -1 if path could not be opened for
    // modification or is not seekable.
    int apply(word_t koop, int ref, string path) {
        // open path for reading and writing
        fstream mod(path);
        if (!mod.is_open())
            return -1;

        // set len to the number of bits in the file
        uintmax_t len;
        {
            mod.seekg(0);
            auto beg = mod.tellg();
            mod.seekg(0, ios::end);
            auto end = mod.tellg();
            if (beg == -1 || end == -1 || end < beg)
                return -1;
            uintmax_t diff = end - beg;         // number of bytes
            len = diff << 3;                    // number of bits
            if ((len >> 3) != diff)
                len = -1;                       // max out on overflow
        }

        // find the longest codeword that will fit (the longest codeword will
        // change the fewest bits)
        if (cat.find(koop) == cat.end())
            cat.insert(make_pair(koop, codewords(koop)));
        auto const& codes = cat.at(koop).codes;
        auto code = codes.upper_bound(len);
        if (code != codes.begin())
            code--;
        if (code->first > len)
            return 0;

        // generate a random bit offset at which to start applying the codeword
        random_device ran;
        default_random_engine gen(ran());
        uniform_int_distribution<uintmax_t> dist(0, len - code->first);
        auto off = dist(gen);

        // exclusive-or the codeword with path starting at bit offset off
        auto const& ones = code->second.ones;
        char buf[1];
        auto cur = (off + ones[0]) >> 3;
        mod.seekg(cur);
        mod.read(buf, 1);
        for (auto one : ones) {
            auto pos = off + one;
            unsigned low = pos & 7;
            pos >>= 3;
            if (cur != pos) {
                mod.seekp(cur);
                mod.write(buf, 1);
                cur = pos;
                mod.seekg(cur);
                mod.read(buf, 1);
            }
            buf[0] ^= ref ? 1 << low : 0x80 >> low;
        }
        mod.seekp(cur);
        mod.write(buf, 1);
        return ones.size();
    }
};

// File to load codewords from.
#define CODES "codewords.txt"

// See the comments at the top for the command line contents.
int main(int argc, char **argv) {
    // load all of Koopman's codewords
    catalog cat;
    int ret = cat.load(CODES);
    if (ret == -1)
        cerr << "could not open " CODES "\n";
    else if (ret)
        cerr << ret << " invalid or too large codeword description" <<
                (ret == 1 ? "" : "s") << " in " CODES "\n";

    // if no arguments, show the number of CRCs and codewords loaded
    if (argc == 1) {
        auto got = cat.counts();
        cout << "loaded " << got.first << " CRC polynomials and " <<
                got.second << " codewords\n";
        return 0;
    }

    // validate the command-line arguments, and convert the provided polynomial
    // into the Koopman convention
    if (argc != 5 ||
        !regex_match(argv[1], regex("\\d+")) ||
        !regex_match(argv[2], regex("[0-9a-fA-F]+")) ||
        !regex_match(argv[3], regex("0|1"))) {
        cerr << "usage: ruse deg poly ref path\n";
        return 1;
    }
    word_t koop;
    try {
        auto n = stoint<unsigned>(argv[1]) - 1;
        koop = stoint<word_t>(argv[2], 0, 16);
        if (n >= numeric_limits<word_t>::digits || (koop & 1) != 1)
            throw out_of_range("");
        koop = ((word_t)1 << n) + (koop >> 1);
    }
    catch (out_of_range const& e) {
        cerr << "arguments out of range\n";
        return 1;
    }

    // apply a selected codeword to the file
    ret = cat.apply(koop, argv[3][0] - '0', argv[4]);
    if (ret == -1)
        cerr << "could not open " << argv[4] << " for modification\n";
    else if (ret == 0)
        cerr << argv[4] << " shorter than the shortest codeword\n";
    else
        cout << "flipped " << ret << " bits in " << argv[4] << "\n";
    return ret <= 0;
}
