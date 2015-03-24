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
        i += encode(c) - '0';
    }

    return i;
}
