#include <sm_common.hpp>

unsigned long int strtob4(const char *str)
{
    register unsigned long int i;
    register const char *s;
    register char c;

    i = 0;
    s = str;

    for (c = *s; c != '\0'; c = *++s) {
        i *= 4;
        i += code[c] - '0';
    }

    return i;
}

// In-place conversion of kmer to its reverse complement.
void krevcomp(char kmer[])
{
    int c, i, j;
    for (i = 0, j = KMER_LEN - 1; i < j; i++, j--)
    {
        c = kmer[i];
        kmer[i] = comp[kmer[j]];
        kmer[j] = comp[c];
    }
}
