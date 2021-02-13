/* flip.c -- use the output of spoof to flip bits in a file

  Copyright (C) 2021 Mark Adler

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
   Apply the output of spoof, read from stdin, to the file whose path is the
   first (and only) command-line argument.

   Usage:

       spoof < request | flip file

   where request is the input to spoof (see spoof.c), and file is the file in
   which to flip the bits chosen by spoof.
 */

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/errno.h>

// Flip bit bit at position pos in file. This defers the read and write to the
// file until pos changes. Call with pos == -1 to do the last exclusive-or. bit
// will be ignored in that case. Return 0 on success, -1 if pos is past the end
// of the file, or a positive errno return value due to a file operation error.
static int flip(FILE *file, off_t pos, int bit) {
    // Keep track of position and exclusive-ors until position changes.
    static int xor = 0;
    static off_t curr = -1;

    // If the position changed, first execute the previous exclusive-or.
    if (curr != pos && curr >= 0) {
        int ret = fseeko(file, curr, SEEK_SET);
        if (ret)
            return errno;
        int got = fgetc(file);
        if (got == EOF)
            return -1;
        got ^= xor;
        ret = fseeko(file, curr, SEEK_SET);
        if (ret)
            return errno;
        ret = fputc(got, file);
        if (ret == EOF)
            return errno;
        xor = 0;
    }
    curr = pos;

    // If this is just a flush, then done.
    if (pos < 0)
        return 0;

    // Accumulate exclusive-or bits.
    xor |= 1 << bit;
    return 0;
}

// Apply the bit flips read from stdin to file. Non-digit characters are
// ignored, except they are considered to terminate a decimal number. Return 0
// on success, -2 if a bit offset is out of range, or a flip() return value.
static int apply(FILE *file) {
    int state = 0, bit, ch;
    off_t pos;
    do {
        ch = getchar();
        int dig = ch >= '0' && ch <= '9' ? ch - '0' : -1;
        switch (state) {
        case 0:
            // looking for start of file position
            if (dig != -1) {
                pos = dig;
                state = 1;
            }
            break;
        case 1:
            // in file position
            if (dig != -1)
                pos = 10 * pos + dig;
            else
                state = 2;
            break;
        case 2:
            // looking for start of bit offset
            if (dig != -1) {
                bit = dig;
                if (bit > 7)
                    return -2;
                state = 3;
            }
            break;
        case 3:
            // in bit offset
            if (dig != -1) {
                bit = 10 * bit + dig;
                if (bit > 7)
                    return -2;
            }
            else {
                // flip bit bit at pos
                int ret = flip(file, pos, bit);
                if (ret)
                    return ret;
                state = 0;
            }
        }
    } while (ch != EOF);

    // flush any remaining exclusive-ors
    return flip(file, -1, 0);
}

// Update the file named on the command line with the output of spoof read from
// stdin.
int main(int argc, char **argv) {
    if (argc != 2) {
        fputs("* expecting one argument: the name of the file\n", stderr);
        return 1;
    }
    FILE *file = fopen(argv[1], "r+b");
    if (file == NULL) {
        fprintf(stderr, "* could not open %s for reading and writing\n",
                argv[1]);
        return 1;
    }
    int ret = apply(file);
    fclose(file);
    if (ret) {
        if (ret == -2)
            fputs("* bit offset out of range -- aborted\n", stderr);
        else if (ret == -1)
            fputs("* file position past end -- aborted\n", stderr);
        else
            fprintf(stderr, "* i/o error %s -- aborted\n", strerror(ret));
        return 1;
    }
    return 0;
}
