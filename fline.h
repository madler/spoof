/* fline -- read lines from a file

  Copyright (C) 2016 Mark Adler

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

// Read delimited lines from a file, fast. This is more than twice as fast as
// BSD fgetln() or POSIX getline(), as measured on Darwin. fline is written to
// be compatible with the C99 standard, not relying on POSIX or other standards
// outside of C99.

#include <stdio.h>      // FILE, size_t

// Type for an fline state. A pointer to this state is returned or used by all
// of the routines here.
typedef struct fline_s fline_t;

// Return a new fline state, a pointer to the fline_t type. This allocates
// memory for the state and assigns the file 'in' to the state, which must be
// readable. Reading will start at in's current file pointer. fline_start()
// returns NULL if it was not able to allocate memory.
fline_t *fline_start(FILE *in);

// Return a line from state, where lines are separated by '\n' characters. The
// returned line includes the '\n', unless it is the last line in the file, and
// the file doesn't end with a '\n'. fline() returns a pointer to the line, and
// sets *len to the length of the line. If fline() returns non-NULL and *len is
// zero, then the end of the file has been reached, or there was a read error.
// ferror(in) can be used to distinguish between those, as for fread(). After
// that, subsequent calls with state will continue to return zero in *len. The
// file provided to the state is not closed.
//
// The returned line is no longer available when fline() or fline_delim() is
// called with the same state to get the next line, or fline_end() is called to
// release the state.
//
// The line is not a C string, i.e. it is not terminated with a null ('\0')
// character, and may in fact have embedded nulls from the file.
//
// The contents of the line may be modified by the caller, within the bounds of
// what was returned. In the cases that the returned line does not have a
// terminating '\n' or a zero length is returned, the caller may modify the
// byte after the returned line. This allows the caller to always either
// replace the terminating '\n' with a null ('\0'), or append a terminating
// null to the last line without a '\n'. The caller should then also consider
// replacing any embedded nulls in the line with some other byte value, in
// order to assure correct use of the line as a C string.
//
// in's file pointer must not be changed between fline() calls.
//
// If fline() returns NULL, then there was a memory allocation error or a
// size_t overflow. In that case, *len is not set, and the state remains valid.
// If it was a memory allocation error, then the operation can be retried after
// more memory is made available. The largest line length that can be assured
// to be readable by fline() without a size_t overflow is (SIZE_MAX >> 2) + 1.
// If size_t is 32 bits, then the assured limit is a 1 GiB line.
char *fline(fline_t *state, size_t *len);

// fline_delim() does what fline() does, but a delimiter other than '\n' can be
// provided as delim. Calls to fline() and fline_delim(), as well as calls to
// fline_delim() with different delimiters can be mixed, simply changing the
// delimiter searched for on each call.
//
// If delim is not in the range 0..255, then the remainder of the entire file
// is returned. For that case, note the returned line length limit in the
// fline() description above.
char *fline_delim(fline_t *state, size_t *len, int delim);

// Return a pointer to and set *len to the length of the unused contents of the
// buffer in state. This can be used to recover unused data read from the file,
// when done reading lines while still in the middle of the file. The returned
// data remains available under the same conditions as for a returned line.
// This does not change the state.
char *fline_remains(fline_t *state, size_t *len);

// Reuse an fline state, assigning the file 'in'. This clears the state, but
// retains the allocated buffer memory (sans contents) from the previous use.
void fline_reuse(fline_t *state, FILE *in);

// End an fline state, freeing the allocated buffer and the state itself. The
// file that was provided for state is not closed.
void fline_end(fline_t *state);
