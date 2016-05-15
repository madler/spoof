// Copyright (C) 2016 Mark Adler
// See fline.h for a description of these functions and the license.

#include <stdio.h>      // size_t, FILE, NULL, fread()
#include <stdlib.h>     // malloc(), free(), realloc()
#include <string.h>     // memchr(), memcpy()
#include "fline.h"

// Internal structure for line reading. This structure is not visible outside
// of this source file.
struct fline_s {
    char *buf;          // input buffer -- returned line is entirely in here
    size_t size;        // allocated size of the input buffer, which can grow
    size_t have;        // index right after end of available data in buf
    size_t pos;         // index of start of available data in buf
    FILE *in;           // input file, set to NULL when EOF is reached
};

void fline_reuse(fline_t *state, FILE *in) {
    state->have = 0;
    state->pos = 0;
    state->in = in;
}

fline_t *fline_start(FILE *in) {
    fline_t *state = malloc(sizeof(fline_t));
    if (state == NULL)
        return NULL;                    // out of memory
    state->size = 32768U;               // must be > 1 and a power of 2
    state->buf = malloc(state->size);
    if (state->buf == NULL) {
        free(state);
        return NULL;                    // out of memory
    }
    fline_reuse(state, in);
    return state;
}

void fline_end(fline_t *state) {
    if (state != NULL) {
        if (state->buf != NULL)
            free(state->buf);
        free(state);
    }
}

char *fline_delim(fline_t *state, size_t *len, int delim) {
    // if the buffer is empty, fill it with input data
    if (state->pos == state->have) {
        if (state->in == NULL) {
            *len = 0;                   // end of file, return 0
            return state->buf;
        }
        state->pos = 0;
        state->have = fread(state->buf, 1, state->size, state->in);
        if (state->have < state->size)
            state->in = NULL;           // end of file reached
    }

    // scan and read input until a delimiter is found or end of file
    size_t next = state->pos;           // start of unscanned input
    for (;;) {
        // scan for the delimiter in the available and unscanned input
        char *hit = delim < 0 || delim > 255 ? NULL :
                    memchr(state->buf + next, delim, state->have - next);

        // if found or if at end of file, return the line
        if (hit != NULL || state->in == NULL) {
            size_t beg = state->pos;
            state->pos = hit == NULL ? state->have : (hit - state->buf) + 1;
            *len = state->pos - beg;
            return state->buf + beg;
        }

        // if the first half of the buffer is now unused, then move the second
        // half of the buffer to the first half (move only the available input)
        next = state->size >> 1;
        if (state->pos >= next) {
            memcpy(state->buf + state->pos - next, state->buf + state->pos,
                   state->size - state->pos);
            state->pos -= next;
        }

        // otherwise allocate more memory, doubling the size of the buffer
        else {
            next = state->size;
            if ((next << 1) == 0)
                return NULL;            // size_t overflow
            char *mem = realloc(state->buf, next << 1);
            if (mem == NULL)
                return NULL;            // out of memory
            state->buf = mem;
            state->size = next << 1;
        }

        // either way, load more input into the second half of the buffer,
        // then search again for a delimiter
        state->have = next + fread(state->buf + next, 1, next, state->in);
        if (state->have < state->size)
            state->in = NULL;           // end of file reached
    }
}

char *fline(fline_t *state, size_t *len) {
    return fline_delim(state, len, '\n');
}

char *fline_remains(fline_t *state, size_t *len) {
    *len = state->have - state->pos;
    return state->buf + state->pos;
}
