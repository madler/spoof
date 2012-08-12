/* spoof.c -- modify a message to have a desired CRC
  version 1.0, August 12th, 2012

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
 */

/*
   Given an n-bit CRC polynomial and n bit locations in a message of specified
   length, determine what to set those bit locations to in order to get a
   specified CRC value -- not all such sets of bit locations have a solution.

   spoof is used by taking a sequence and its n-bit CRC value, selecting n bit
   locations in the sequence to potentially change, and exclusive-oring the CRC
   value with the desired CRC.  The bit locations and that difference between
   the CRCs is provided to spoof.  Then spoof delivers a subset of the bit
   locations that are to be inverted (0 -> 1 or 1 -> 0).  Upon inverting, that
   sequence now has the desired CRC.  If spoof reports that that set of bit
   locations has no solution, then a different set of bit locations can be
   tried by the user.

   The input is read from stdin.  The format of the input is:

     dimension reflect polynomial
     length
     offset_1 position_1
     offset_2 position_2
     ...
     offset_n position_n
     crc

   The first line describes the CRC, where 'dimension' is the number of bits in
   the crc in decimal, 'reflect' is 1 for a reflected crc or 0 for a
   non-reflected crc, and 'polynomial' is the crc polynomial in hexadecimal.
   The polynomial is represented by its low coefficients (i.e. not including
   the coefficent of x^dimension, which is always 1), with the x^0 coefficient
   placed in the least significant bit for a non-reflected CRC, or in the most
   significant bit (of a dimension-bits word) for a reflected CRC.  Reflection
   of the CRC is applied on both input and output.  There is no specification
   required for pre or post processing of the CRC, since the result of spoof is
   independent of such processing.

   On the next line 'length' is the length of the sequence in bytes, expressed
   in decimal, where each byte is eight bits.  Then there are n bit locations,
   where n is equal to dimension.  Each bit location consists of 'offset',
   which is the distance of the location in bytes from the start of the
   sequence, in decimal, where zero is the first byte in the sequence, and
   'position' which is the location of the bit in the byte in decimal, with
   zero representing the least-significant bit.  Lastly 'crc' is the
   exclusive-or of the initial and desired CRCs, expressed in hexadecimal.  Any
   blank character can be used to separate the values.  New line characters can
   be used as shown above for readability, but are not required.

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
     89
     37 0
     41 0
     45 0
     49 0
     f

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
   changed bits are not consequential to the transmitted message.

   We will place in each candidate bit location in D a variable, named x_0,
   x_1, etc., with all of the other bits in D set to zero.  The equation: r =
   crc(D) for an n-bit CRC can be seen as n binary equations in the x_i, over
   GF(2).  We will define n such locations x_i, since then we have n equations
   in n unknowns, which is a minimal set of candidate locations which will have
   at most one solution.

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
   inverted to get a sequence B that has the CRC q.  If M is singular, there is
   no solution for the given set of bit locations (regardless of r).  The user
   can then try a different set of bit locations.

   The described application of spoof works as well for CRC's calculated with
   pre and/or post-processing, where the initial CRC value may be non-zero, and
   the final CRC value may be exclusive-or'ed with a constant.  That processing
   can be seen as simply exclusive-or'ing a single constant with the CRC, where
   that constant depends only on the length.  spoof does its calculations using
   only a "pure" CRC with no pre- or post-processing.  When the bits specified
   by spoof are inverted, the same constant is still part of the original CRC,
   and so its effect is unchanged in the resulting CRC.

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

/* CRC description (with no pre or post processing) */
typedef struct {
    unsigned short dim; /* number of bits in CRC */
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
local inline void *alloc(size_t size)
{
    void *space;

    space = malloc(size);
    if (space == NULL)
        fail("out of memory");
    return space;
}

/* Run the low eight bits in val through crc model. */
local inline word crc_byte(word crc, unsigned val, model_t model)
{
    unsigned k;

    if (model.ref) {
        crc ^= val & 0xff;
        for (k = 0; k < 8; k++)
            crc = crc & 1 ? (crc >> 1) ^ model.poly : crc >> 1;
    }
    else {
        word m;

        m = (word)1 << (model.dim - 1);
        crc ^= ((word)val & 0xff) << (model.dim - 8);
        for (k = 0; k < 8; k++)
            crc = crc & m ? (crc << 1) ^ model.poly : crc << 1;
        crc &= ((word)1 << model.dim) - 1;
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
local void gf2_matrix_square(word *square, const word *mat, unsigned dim)
{
    unsigned n;

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
    static word *power[sizeof(range) * 8];

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
            unsigned n;
            word row;

            /* check and set state, allocate space for first two operators */
            first = model;
            power[0] = alloc(model.dim * sizeof(word));
            power[1] = alloc(model.dim * sizeof(word));

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
        power[have] = alloc(model.dim * sizeof(word));
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

/* Compute the inverse of mat, with the result inv.  Return 0 on success, 1 if
   mat is singular.  dim is the dimension of the matrices. */
local int gf2_matrix_invert(word *inv, const word *mat, unsigned dim)
{
    unsigned k, j;
    word pos, a[dim];

    /* copy mat to local storage and create adjoining identity matrix */
    for (k = 0, pos = 1; k < dim; k++, pos <<= 1) {
        a[k] = mat[k];
        inv[k] = pos;
    }

    /* make a[] the identity matrix using column swaps and column subtractions
       (exclusive-or), and perform the same operations on inv[] -- then inv[]
       will be the inverse */
    for (k = 0, pos = 1; k < dim; k++, pos <<= 1) {
        /* find a subsequent row where column k is 1, make that row k with a
           swap if necessary -- if there isn't any such row, then the starting
           matrix is singular, in which case return an error */
        if ((a[k] & pos) == 0) {
            word tmp;

            for (j = k + 1; j < dim; j++)
                if (a[j] & pos)
                    break;
            if (j == dim)           /* no such row, matrix is singular */
                return 1;
            tmp = a[k], a[k] = a[j], a[j] = tmp;
            tmp = inv[k], inv[k] = inv[j], inv[j] = tmp;
        }

        /* subtract row k from all the other rows with a 1 in that column */
        for (j = 0; j < dim; j++)
            if (j != k && (a[j] & pos) != 0) {
                a[j] ^= a[k];
                inv[j] ^= inv[k];
            }
    }
    return 0;
}

/* Compute the crc of a sparse sequence with 1's at loci[0..num-1] (assumed to
   be sorted by offset).  Any that have pos > 7 are ignored. */
local word crc_sparse(const struct locus *loci, unsigned num, range len,
                      model_t model)
{
    unsigned k;         /* index of loci */
    unsigned val = 0;   /* sequence byte consisting of one or more ones */
    word crc = 0;       /* computed crc */
    range at = 0;       /* crc calculation at this offset so far */

    /* go through each location, deferring the use of val in case a byte has
       more than one location in it */
    for (k = 0; k < num; k++)
        if (loci[k].pos < 8) {                  /* skip marked locations */
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

/* Solve for the set of loci and the desired crc.  Return true if there is no
   solution, otherwise return the solution in sol.  The solution is a word with
   1's in the positions corresponding to loci that need to be inverted. */
local int crc_solve(word *sol, const struct locus *loci, range len, word want,
                    model_t model)
{
    unsigned k;
    word mat[model.dim], inv[model.dim];

    /* for each bit position, calculate the crc of the sequence of len zero
       bytes except for a single 1 bit at that bit position */
    for (k = 0; k < model.dim; k++)
        mat[k] = crc_sparse(loci + k, 1, len, model);

    /* invert the resulting matrix of crc's (return if singular) */
    k = gf2_matrix_invert(inv, mat, model.dim);
    if (k)
        return 1;

    /* apply inverse to desired crc to get solution, which indicates which bits
       in loci to set */
    *sol = gf2_matrix_times(inv, want);
    return 0;
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
    int k;

    k = 0;
    do {
        n /= 10;
        k++;
    } while (n);
    return k;
}

/* Read sequence length, bit positions, and desired crc difference from stdin.
   Compute and display the solution, which is a subset of the provided bit
   positions to invert in the sequence. */
int main(void)
{
    int ret;                    /* general function return value */
    unsigned k;                 /* counter for locations, bits */
    FILE *in = stdin;           /* input file */
    model_t model;              /* CRC model */
    range len;                  /* sequence length */
    struct locus *loci;         /* bit locations */
    word want;                  /* desired crc */
    word sol;                   /* solution for bit locations */
    word crc;                   /* calculated crc to check solution */
    char *just;                 /* string of spaces for justification */

    /* read and validate the input file */
    for (k = 0, crc = 1; crc; k++, crc <<= 1)
        ;
    ret = fscanf(in, " %hu %hd %" WORDFMT,   /* crc description */
                 &model.dim, &model.ref, &model.poly);
    if (ret == 3 && model.dim > k)
        fail("CRC too long for crc integer type spoof was compiled with");
    if (ret < 3 || model.dim == 0 || model.ref < 0 || model.ref > 1 ||
                   model.poly > ((word)1 << model.dim) - 1)
        fail("invalid CRC description");
    if ((model.poly & ((word)1 << (model.ref ? model.dim - 1 : 0))) == 0)
        fail("invalid polynomial (you may need to reverse the bits)");
    loci = alloc(model.dim * sizeof(struct locus));
    ret = fscanf(in, " %" RANGEFMT, &len);      /* length of sequence */
    if (ret < 1 || len < ((model.dim + 7) >> 3))
        fail("invalid sequence length (must be at least length of CRC)");
    for (k = 0; k < model.dim; k++) {
        ret = fscanf(in, " %" RANGEFMT " %hd",   /* offset and bit position */
                     &(loci[k].off), &(loci[k].pos));
        if (ret < 2 || loci[k].off >= len || loci[k].pos > 7)
            fail("invalid bit location specification");
    }
    ret = fscanf(in, " %" WORDFMT, &want);      /* desired crc difference */
    if (ret < 1 || want > ((word)1 << model.dim) - 1)
        fail("invalid target CRC");
    while ((ret = getc(in)) != EOF)
        if (ret != ' ' && ret != '\t' && ret != '\n' && ret != '\r' &&
            ret != '\v' && ret != '\f') {
            warn("junk at end of input ignored");
            break;
        }

    /* sort the bit locations for the benefit of crc_sparse() -- not needed
       for crc_solve(), but it doesn't hurt to have a sorted result */
    mergesort(loci, (int)(model.dim), sizeof(struct locus), locus_order);

    /* solve for the values of the given bit locations to get want */
    ret = crc_solve(&sol, loci, len, want, model);
    if (ret)
        fail("no solution for the given bit locations -- try another set");

    /* invalidate the positions of the zeros to remove them from the check and
       from the output */
    for (k = 0, crc = 1; k < model.dim; k++, crc <<= 1)
        if ((sol & crc) == 0)
            loci[k].pos += 8;

    /* check the crc of a sequence with ones at the remaining locations */
    crc = crc_sparse(loci, model.dim, len, model);
    if (want != crc)
        fail("internal algorithm error");

    /* determine the number of digits in the largest offset and create a string
       with leading blanks for formatting the output nicely */
    for (k = model.dim; k > 0; k--)
        if (loci[k - 1].pos < 8)
            break;
    if (k == 0)
        fail("internal algorithm error");
    ret = decimal_digits(loci[k - 1].off);
    if (ret < 6)
        ret = 6;
    just = alloc(ret - 5);
    memset(just, ' ', ret - 6);
    just[ret - 6] = 0;

    /* output what bits to invert to get the desired crc */
    puts("invert these bits in the sequence:");
    printf("%s%s\n", just, "offset bit");
    free(just);
    for (k = 0; k < model.dim; k++)
        if (loci[k].pos < 8)
            printf("%*" RANGEFMT " %d\n", ret, loci[k].off, loci[k].pos);

    /* clean up and return success */
    crc_zeros_operator(-1, model);
    free(loci);
    return 0;
}
