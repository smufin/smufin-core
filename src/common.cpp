/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. You
 * should have received a copy of the SMUFIN Public License along with this
 * file. If not, see <https://github.com/smufin/smufin-core/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 *
 */

#include "common.hpp"

// S tring to integer conversion. Parses a null-terminated string, interpreting
// its content as an integral number in base 4.
uint64_t strtob4(const char *str)
{
    register uint64_t i = 0;
    register const char *s = str;
    register char c;
    for (c = *s; c != '\0'; c = *++s) {
        i *= 4;
        i += sm::code[c] - '0';
    }
    return i;
}

// Coded base 4 sequence to string conversion.
void b4tostr(uint64_t code, int len, char *str)
{
    for (int i = 0; i < len; i++) {
        str[len - i - 1] = sm::alpha[code & 3];
        code >>= 2;
    }
    str[len] = '\0';
}

// Low-quality phred score counter. Returns number of bases with a quality
// score below 20.
int lq_count(const char *str, int len)
{
    int lq = 0;
    for (int i = 0; i < len; i++) {
        int phred = str[i] - 33;
        if (phred < 20)
            lq++;
    }
    return lq;
}

// In-place conversion of a sequence to its reverse.
void rev(char seq[], int len)
{
    int i, j;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        const char c = seq[i];
        seq[i] = seq[j];
        seq[j] = c;
    }
}

// In-place conversion of a sequence to its reverse complement.
void revcomp(char seq[], int len)
{
    int i, j;
    for (i = 0, j = len - 1; i <= j; i++, j--) {
        const char c = sm::comp[seq[i]];
        seq[i] = sm::comp[seq[j]];
        seq[j] = c;
    }
}
