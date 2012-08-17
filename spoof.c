/* spoof.c -- modify a message to have a desired CRC
  version 1.1, August 16th, 2012

  Copyright (C) 2012 Mark Adler

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

/* Version History:
   1.0    12 Aug 2012  First version
   1.1    16 Aug 2012  Allow underconstrained set to avoid singular matrices
 */

/*
   Given a k-bit CRC polynomial and n >= k bit locations in a message of
   specified length, determine what to set those bit locations to in order to
   get a specified CRC value.  Not all such sets of bit locations have a
   solution, but providing n > k candidate bit locations to change reduces the
   probability of no solution.

   spoof is used by taking a sequence and its CRC value, selecting n bit
   locations in the sequence to potentially change, and exclusive-oring the CRC
   value with the desired CRC.  The bit locations and that difference between
   the CRCs is provided to spoof.  Then spoof delivers a subset of the bit
   locations that are to be inverted (0 -> 1 or 1 -> 0).  Upon inverting, that
   sequence now has the desired CRC.  If spoof reports that that set of bit
   locations has no solution, then a different or larger set of bit locations
   can be tried by the user.

   The input is read from stdin.  The format of the input is:

     dimension reflect polynomial
     crc
     length
     offset_1 position_1
     offset_2 position_2
     ...
     offset_n position_n

   The first line describes the CRC, where 'dimension' is the number of bits in
   the crc in decimal, 'reflect' is 1 for a reflected crc or 0 for a non-
   reflected crc, and 'polynomial' is the crc polynomial in hexadecimal. The
   polynomial is represented by its low coefficients (i.e. not including the
   coefficent of x^dimension, which is always 1), with the x^0 coefficient
   placed in the least significant bit for a non-reflected CRC, or in the most
   significant bit (of a dimension-bits word) for a reflected CRC.  Reflection
   of the CRC is applied on both input and output.  There is no specification
   required for pre or post processing of the CRC, since the result of spoof is
   independent of such processing.

   On the next line 'crc' is the exclusive-or of the initial and desired CRCs,
   expressed in hexadecimal. 'length' is the length of the sequence in bytes,
   expressed in decimal, where each byte is eight bits.  Then there are n bit
   locations, where n is equal to or greater than dimension.  Each bit location
   consists of 'offset', which is the distance of the location in bytes from
   the start of the sequence, in decimal, where zero is the first byte in the
   sequence, and 'position' which is the location of the bit in the byte in
   decimal, with zero representing the least-significant bit.  'offset' must be
   less than 'length', and 'position' must be less than eight.  The end of the
   list is indicated by the end of the input file.  Any blank character can be
   used to separate the values.  New line characters can be used as shown above
   for readability, but are not required.

   Some examples for <dimension reflect polynomial> for common CRCs are:

     32 1 edb88320          ZIP/GZIP/PNG
     32 0 04c11db7          BZIP2/POSIX/MPEG2 (same polynomial as ZIP)
     16 1 8408              X.25/KERMIT/HDLC/CCITT
     64 1 c96c5795d7870f42  XZ

   If the sequence of message bits is not a multiple of eight, prepend the
   sequence with zero bits until it is, and don't specify any locations in the
   prepended bits.  Then compensate for the number of prepended bits when
   interpreting the output of spoof.

   The output of spoof is written to stdout in readable form, as a table of
   offset and position pairs, one pair per line, that should be inverted in the
   sequence, preceded by two lines of instruction and table header.  E.g.:

     invert these bits in the sequence:
     offset bit
         33 1
         36 2

   These pairs will be a subset of the pairs provided in the input.

   An example of a complete input file, using CRC-4/ITU (a four-bit CRC) is:

     4 1 c
     f
     89
     37 0
     41 0
     45 0
     49 0

   The resulting output is:

     invert these bits in the sequence:
     offset bit
         41 0

   The execution time of spoof is proportional to log(length).  So spoof can be
   used for extremely long sequences and still return a solution very rapidly.
 */

/*
   How it works:

   Given two sequences of the same length, the CRC of the exclusive-or of the
   two sequences is equal to the exclusive-or of the CRCs of the sequences
   separately.  This relation is a consequence of the linearity of the CRC over
   the Galois field of order two, referred to as GF(2).  GF(2) consists of just
   the two elements 0 and 1, and the operations exclusive-or and logical-and,
   which take the place of arithmetic's addition and multiplication operations
   respectively.  This additive or superposition property allows spoof to never
   need to know the message contents in order to find a solution.  All it needs
   is the before and after CRCs, or really just the exclusive-or of those two
   CRCs.

   Given a sequence A and CRC p, we would like to modify A to a new sequence B,
   to give a specified CRC q.  So {A, p} -> {B, q}.  We are given A, p, and q,
   and we need to find B.

   There are many answers for B.  In order to narrow those down, we would like
   to make only a small number of changes to A.  Let D = A ^ B and r = p ^ q.
   We have from the above that r = crc(D).  We would like for D to be mostly
   zeros, with just a small number of ones, which represent the number of bit
   locations where A and B differ.  r is simply calculated from p and q, which
   are known.  We will pick a set of bit locations in D that we will allow
   spoof to modify.  These bit locations can be anywhere, such as all grouped
   at the end or beginning, randomly scattered in the sequence, the low bits of
   selected insignificant decimal digits, or perhaps other choices where the
   changed bits are not consequential to the transmitted message.  spoof can
   also be used to attempt to correct a set of known erasure locations using
   the CRC.

   We will place in each candidate bit location in D a variable, named x_0,
   x_1, etc., with all of the other bits in D set to zero.  The equation: r =
   crc(D) for a k-bit CRC can be seen as k binary equations in the x_i, over
   GF(2).  We will define n such locations x_i, where n >= k, since then we
   have at k equations with at least k unknowns.  Out of the n x_i, we will
   look for a subset k x_i that results in a solution.

   Given the length of the sequence, r, and the locations of the x_i, spoof
   will determine the values of the x_i, from which D can be constructed.  Then
   B = A ^ D, where q = crc(B).  Or more simply, for each x_i that is one,
   invert the bit at that location in A to get B.  spoof does not need to know
   A, just the locations of the x_i.

   For each x_i, we consider a sequence X_i which is all zeros except for a
   single one at the x_i location.  We then calculate the CRC of each X_i,
   giving c_i = CRC(X_i).  We now have n c_i values.  If there is a solution,
   then there is a subset of the c_i that, when exclusive-ored together, is
   equal to r.  To solve, we construct the matrix M that consists of the
   columns c_i.  If x is the vector x_i, then we have M x = r.  We take the
   inverse of M, which if it exists, gives the solution x = Inverse(M) r.  For
   the x_i with the value one, the corresponding locations in A need to be
   inverted to get a sequence B that has the CRC q.  If all square subset of M
   are singular, there is no solution for the given set of bit locations
   (regardless of r).  The user can then try a different or larger set of bit
   locations.

   The described application of spoof works as well for CRC's calculated with
   pre and/or post-processing, where the initial CRC value may be non-zero, and
   the final CRC value may be exclusive-or'ed with a constant.  That processing
   can be seen as simply exclusive-or'ing a single constant with the CRC, where
   that constant depends only on the length of the sequence the CRC is over.
   spoof does its calculations using only a "pure" CRC with no pre- or
   post-processing.  This is permitted since spoof is provided the exlusive-or
   of two sequences of the same length, which cancels exclusive-or'ed constant,
   leaving the pure CRC of the two sequences exclusive-or'ed.

   The usual way to calculate c_i = crc(X_i) takes an amount of time linear in
   the length of the sequence.  However for sparse sequences, that execution
   time can be shortened dramatically by constructing matrix operators that
   represent the application of a series of zeros to the CRC.  We construct a
   matrix representing the effect on the CRC of running a single zero bit
   through the CRC.  Call it Z.  Then we successively square that matrix to get
   operators for more zeros.  Z**8 represents running a byte of zeros thruogh
   the CRC (where ** means to the power of).  Z**16 is two bytes.  Z**32 is
   four bytes.  And so on.  Then we simply decompose the length into a sum of
   powers of two, and apply the corresponding operators for those numbers of
   zeros to the CRC.

   As a result, spoof runs in O(log n) time, where n is the length of the
   sequence being spoofed.  The execution time also depends on the dimension of
   the CRC in order to square matrices.  Let d be that dimension.  Then spoof
   runs in O(d**2 log(n)) time.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define local static

/* Types to use for CRC's and sequence lengths and offsets.  In general these
   should be the largest integer types available to maximize the problems that
   can be solved.  word could be made a smaller type if speed is paramount and
   the size of the word type is known to cover the CRC polynomials that will be
   presented.
 */
typedef unsigned long long word;    /* unsigned type for crc values */
typedef unsigned long long range;   /* unsigned type for sequence offsets */
#define WORDFMT "llx"               /* printf / scanf format for word (hex) */
#define RANGEFMT "llu"              /* printf / scanf format for range */
#define WORDBITS ((int)sizeof(word)<<3)
#define ONES(n) ((n) >= WORDBITS ? (word)0 - 1 : ((word)1 << (n)) - 1)

/* CRC description (with no pre or post processing) */
typedef struct {
    short dim;          /* number of bits in CRC */
    short ref;          /* if true, bit-reflected input and output */
    word poly;          /* polynomial representation (ordered per ref) */
} model_t;

/* Location of a bit that can be modified to get the desired CRC. */
struct locus {
    range off;          /* byte offset in sequence */
    short pos;          /* position in byte (0..7) */
};

/* Issue error message (all error messages go through here). */
local inline void warn(const char *why)
{
    fprintf(stderr, "spoof: %s\n", why);
}

/* Fail and exit with error message. */
local inline void fail(const char *why)
{
    warn(why);
    exit(1);
}

/* Assured memory allocation. */
local inline void *alloc(void *space, size_t size)
{
    space = realloc(space, size);
    if (space == NULL)
        fail("out of memory");
    return space;
}

/* Run the low eight bits in val through a crc using model. */
local word crc_byte(word crc, unsigned val, model_t model)
{
    if (model.ref) {
        if (model.dim < 8) {
            crc = (crc ^ val) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 1) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 2) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 3) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 4) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 5) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 6) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = (crc ^ val >> 7) & 1 ? crc >> 1 ^ model.poly : crc >> 1;
        }
        else {
            crc ^= val & 0xff;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
            crc = crc & 1 ? crc >> 1 ^ model.poly : crc >> 1;
        }
    }
    else {
        word mask = (word)1 << model.dim - 1;
        if (model.dim < 8) {
            val <<= model.dim - 1;
            crc = (crc ^ val >> 7) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 6) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 5) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 4) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 3) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 2) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val >> 1) & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = (crc ^ val) & mask ? crc << 1 ^ model.poly : crc << 1;
        }
        else {
            crc ^= (word)val << model.dim - 8;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
            crc = crc & mask ? crc << 1 ^ model.poly : crc << 1;
        }
        crc &= ONES(model.dim);
    }
    return crc;
}

/* Multiply the GF(2) vector vec by the GF(2) matrix mat, returning the
   resulting vector.  The vector is stored as bits in a word, and the matrix
   is similarly stored as words, where the number of words is at least enough
   to cover the position of the most significant 1 bit in the vector (so a
   dimension parameter is not needed). */
local inline word gf2_matrix_times(const word *mat, word vec)
{
    word sum;

    sum = 0;
    while (vec) {
        if (vec & 1)
            sum ^= *mat;
        vec >>= 1;
        mat++;
    }
    return sum;
}

/* Multiply the matrix mat by itself, returning the result in square. dim is
   the dimension of the matrices, i.e., the number of bits in each word (rows),
   and the number of words (columns). */
local void gf2_matrix_square(word *square, const word *mat, int dim)
{
    int n;

    for (n = 0; n < dim; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}

/* Return a matrix that when multiplied by the starting crc is equivalent to
   running 2^k zero bytes through the crc calculation.  The matrices are
   retained in static and allocated storage, so that they are only calculated
   once.  Call crc_zeros_operator(-1, model) to free the allocated storage and
   clear the table.  This routine is not thread-safe. */
local const word *crc_zeros_operator(int k, model_t model)
{
    static int have = 0;
    static model_t first;
    static word *power[sizeof(range) << 3];

    /* if requested or required, release and clear the operator table */
    if (k < 0 || model.dim != first.dim || model.ref != first.ref ||
                 model.poly != first.poly) {
        while (have)
            free(power[--have]);
        if (k < 0)
            return 0;
    }

    /* if necessary, square up to the requested operator */
    while (k >= have) {
        /* first time in: create first two operators (1 and 2 zero bytes) */
        if (have == 0) {
            int n;
            word row;

            /* check and set state, allocate space for first two operators */
            first = model;
            power[0] = alloc(NULL, model.dim * sizeof(word));
            power[1] = alloc(NULL, model.dim * sizeof(word));

            /* generate operator for one zero bit using crc polynomial */
            if (model.ref) {
                power[1][0] = model.poly;
                for (n = 1, row = 1; n < model.dim; n++, row <<= 1)
                    power[1][n] = row;
            }
            else {
                for (n = 0, row = 2; n < model.dim - 1; n++, row <<= 1)
                    power[1][n] = row;
                power[1][n] = model.poly;
            }

            /* square that until we get the operator for eight zero bits */
            gf2_matrix_square(power[0], power[1], model.dim);
            gf2_matrix_square(power[1], power[0], model.dim);
            gf2_matrix_square(power[0], power[1], model.dim);

            /* since we have already allocated the space for it, compute
               the operator for two zero bytes (16 zero bits) */
            gf2_matrix_square(power[1], power[0], model.dim);
            have = 2;
            continue;
        }

        /* square the highest operator so far and put in allocated space */
        power[have] = alloc(NULL, model.dim * sizeof(word));
        gf2_matrix_square(power[have], power[have - 1], model.dim);
        have++;
    }

    /* return the requested operator */
    return power[k];
}

/* Efficiently apply len zero bytes to crc, returning the resulting crc.  The
   execution time of this routine is proportional to log(len).  model is the
   crc description. */
local word crc_zeros(word crc, range len, model_t model)
{
    int n;

    /* apply len zeros to crc */
    if (crc)
        for (n = 0; len; len >>= 1, n++)
            if (len & 1)
                crc = gf2_matrix_times(crc_zeros_operator(n, model), crc);
    return crc;
}

/* Compute the crc of a sparse sequence with 1's at loci[0..num-1] (assumed to
   be sorted by offset). */
local word crc_sparse(const struct locus *loci, int locs, range len,
                      model_t model)
{
    int k;              /* index of loci */
    unsigned val = 0;   /* sequence byte consisting of one or more ones */
    word crc = 0;       /* computed crc */
    range at = 0;       /* crc calculation at this offset so far */

    /* go through each location, deferring the use of val in case a byte has
       more than one location in it */
    for (k = 0; k < locs; k++) {
        /* if at new location, do crc through val if val has ones */
        if (val && loci[k].off != at) {
            crc = crc_byte(crc, val, model);
            at++;
            val = 0;
        }
        /* run zeros through crc up to current location */
        crc = crc_zeros(crc, loci[k].off - at, model);
        at = loci[k].off;
        val |= 1 << loci[k].pos;            /* add a one bit to val */
    }
    /* take care of leftover bits in val, if any */
    if (val) {
        crc = crc_byte(crc, val, model);
        at++;
    }

    /* take care of leftover zeros to run through, return result */
    return crc_zeros(crc, len - at, model);
}

/* Solve M x = c for x, return 0 on success, 1 on failure (singular).  This
   works for rectangluar M as well (cols > rows), where a subset of the x
   values are selected that result in a non-singular square M' over that
   subset.  rows is limited to the number of bits in the word type.  cols is
   not limited (except by stack space). */
local int gf2_matrix_solve(word *x, const word *M, word c, int rows, int cols)
{
    int n = (cols + WORDBITS - 1) / WORDBITS;   /* words to hold cols bits */
    int k;              /* index through columns */
    int j;              /* index through rows */
    int i;              /* index through n words holding cols bits */
    word pos;           /* word with one bit set for current row or column */
    word a[cols];       /* starting matrix, evolving to identity matrix */
    word inv[cols][n];  /* identity matrix, evolving to inverse matrix */

    /* copy mat to local storage and create adjoining identity matrix */
    for (k = 0, j = 0, pos = 1; k < cols; k++, pos <<= 1) {
        if (pos == 0) {
            pos = 1;
            j++;
        }
        a[k] = M[k];
        for (i = 0; i < n; i++)
            inv[k][i] = i == j ? pos : 0;
    }

    /* make a[] the identity matrix using column swaps and column subtractions
       (exclusive-or), and perform the same operations on inv[] -- then inv[]
       will be the inverse */
    for (j = 0, pos = 1; j < rows; j++, pos <<= 1) {
        /* find a subsequent row where column j is 1, make that row j with a
           swap if necessary -- if there isn't any such row, then there is no
           non-singular subset of M, in which case return an error */
        if ((a[j] & pos) == 0) {
            word tmp;

            for (k = j + 1; k < cols; k++)
                if (a[k] & pos)
                    break;
            if (k == cols)          /* no such row, matrix is singular */
                return 1;
            tmp = a[j], a[j] = a[k], a[k] = tmp;
            for (i = 0; i < n; i++)
                tmp = inv[j][i], inv[j][i] = inv[k][i], inv[k][i] = tmp;
        }

        /* subtract row j from all the other rows with a 1 in that column */
        for (k = 0; k < cols; k++)
            if (k != j && (a[k] & pos) != 0) {
                a[k] ^= a[j];
                for (i = 0; i < n; i++)
                    inv[k][i] ^= inv[j][i];
            }
    }

    /* Multiply inverse by c to get result x */
    assert(c <= ONES(rows));
    for (i = 0; i < n; i++)
        x[i] = 0;
    for (j = 0; c; c >>= 1, j++)
        if (c & 1) {
            for (i = 0; i < n; i++)
                x[i] ^= inv[j][i];
        }
    return 0;
}

/* Solve for the set of loci and the desired crc.  Return the number of
   locations to invert, or -1 if there is no solution.  The locations to invert
   are moved to the beginning of loci.  If there is no solution, loci is
   not modified. */
local int crc_solve(struct locus *loci, int locs, range len, word want,
                    model_t model)
{
    int n, k, i;
    word p, sol[(locs + WORDBITS - 1) / WORDBITS];
    word mat[locs];

    /* protect against improper input causing array overruns */
    assert(locs >= model.dim);
    assert(want <= ONES(model.dim));

    /* for each bit position, calculate the crc of the sequence of len zero
       bytes except for a single 1 bit at that bit position */
    for (k = 0; k < locs; k++)
        mat[k] = crc_sparse(loci + k, 1, len, model);

    /* solve mat . sol = want for sol (return if all square subsets of mat are
       singular) */
    k = gf2_matrix_solve(sol, mat, want, model.dim, locs);
    if (k)
        return -1;

    /* move the locations to invert up to the front of loci */
    for (k = 0, n = 0, i = 0, p = 1; k < locs; k++, p <<= 1) {
        if (p == 0) {
            p = 1;
            i++;
        }
        if (sol[i] & p) {
            loci[n] = loci[k];
            n++;
        }
    }
    return n;
}

/* Comparison function for sorting loci, used by mergesort(). */
local int locus_order(const void *a, const void *b)
{
    const struct locus *p = a, *q = b;

    if (p->off != q->off)
        return p->off < q->off ? -1 : 1;
    return p->pos < q->pos ? -1 : (p->pos > q->pos ? 1 : 0);
}

/* Return the number of decimal digits in the unsigned number n. */
local inline int decimal_digits(range n)
{
    int i;

    i = 0;
    do {
        n /= 10;
        i++;
    } while (n);
    return i;
}

/* Read sequence length, bit positions, and desired crc difference from stdin.
   Compute and display the solution, which is a subset of the provided bit
   positions to invert in the sequence. */
int main(void)
{
    int k;                      /* counter for locations, bits */
    word crc;                   /* calculated crc to check solution */
    int ret;                    /* general function return value */
    FILE *in = stdin;           /* input file */
    model_t model;              /* CRC model */
    word want;                  /* desired crc */
    range len;                  /* sequence length */
    struct locus *loci;         /* bit locations */
    int locs;                   /* number of bit locations to look at */
    int flips;                  /* number of bit locations to invert */
    char *just;                 /* string of spaces for justification */

    /* read and validate the input file */
    for (k = 0, crc = 1; crc; k++, crc <<= 1)
        ;
    ret = fscanf(in, " %hd %hd %" WORDFMT,   /* crc description */
                 &model.dim, &model.ref, &model.poly);
    if (ret == 3 && model.dim > k)
        fail("CRC too long for crc integer type spoof was compiled with");
    if (ret < 3 || model.dim < 1 || model.ref < 0 || model.ref > 1 ||
                   model.poly > ONES(model.dim))
        fail("invalid CRC description");
    if ((model.poly & ((word)1 << (model.ref ? model.dim - 1 : 0))) == 0)
        fail("invalid polynomial (you may need to reverse the bits)");
    ret = fscanf(in, " %" WORDFMT, &want);      /* desired crc difference */
    if (ret < 1 || want > ONES(model.dim))
        fail("invalid target CRC");
    ret = fscanf(in, " %" RANGEFMT, &len);      /* length of sequence */
    if (ret < 1 || len < (range)((model.dim + 7) >> 3))
        fail("invalid sequence length (must be at least length of CRC)");
    k = model.dim << 1;
    loci = alloc(NULL, k * sizeof(struct locus));
    locs = 0;
    for (;;) {
        ret = fscanf(in, " %" RANGEFMT " %hd",   /* offset and bit position */
                     &(loci[locs].off), &(loci[locs].pos));
        if (ret == EOF || ret == 0)
            break;
        if (ret < 2 || loci[locs].off >= len ||
                loci[locs].pos < 0 || loci[locs].pos > 7)
            fail("invalid bit location specification");
        locs++;
        if (locs == k) {
            k <<= 1;
            loci = alloc(loci, k * sizeof(struct locus));
        }
    }
    if (locs < model.dim)
        fail("need at least n bit locations for an n-bit CRC");
    if (ret == 0)
        warn("junk at end of input ignored");
    loci = alloc(loci, locs * sizeof(struct locus));

    /* solve for the values of the given bit locations to get want */
    flips = crc_solve(loci, locs, len, want, model);
    if (flips == -1)
        fail("no solution -- try more or different bit locations");

    /* check the crc of a sequence with ones at the given locations -- sort the
       locations by offset first, since crc_sparse() requires that */
    mergesort(loci, flips, sizeof(struct locus), locus_order);
    crc = crc_sparse(loci, flips, len, model);
    if (want != crc)
        fail("internal algorithm error");

    /* output what bits to invert to get the desired crc */
    if (flips) {
        puts("invert these bits in the sequence:");
        ret = decimal_digits(loci[flips - 1].off);
        if (ret < 6)
            ret = 6;
        just = alloc(NULL, ret - 5);
        memset(just, ' ', ret - 6);
        just[ret - 6] = 0;
        printf("%soffset bit\n", just);
        free(just);
        for (k = 0; k < flips; k++)
            printf("%*" RANGEFMT " %d\n", ret, loci[k].off, loci[k].pos);
    }
    else
        puts("no need to invert any bits in sequence");

    /* clean up and return success */
    crc_zeros_operator(-1, model);
    free(loci);
    return 0;
}
