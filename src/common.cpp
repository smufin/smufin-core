#include "common.hpp"

// String to integer conversion. Parses a null-terminated string, interpreting
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

// In-place conversion of sequence to its reverse complement.
void revcomp(char seq[], int len)
{
    int c, i, j;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        c = seq[i];
        seq[i] = sm::comp[seq[j]];
        seq[j] = sm::comp[c];
    }
}
