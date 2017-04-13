# Performance Tuning and Limitations

*smufin* can be easily adapted to work for different kinds of inputs, read
lengths, coverages, etc. This document discusses tradeoffs and assumptions
made for default configurations, as well as hardcoded limitations that would
require additional changes.

## Length of Kmers

The length of kmers, defined in the configuration file as `core.k` is expected
to be in the `[5,32]` range. At the lower end, a minimum of 5 is required
because that's the number of bases that are currently used for partitioning,
although 20 or more is recommended to increase uniqueness. At the higher end,
kmers are limited to 32 mostly due to memory constraints during the `count`
stage. Sequences of up to 32 characters in the 4-base ACGT alphabet can be
stored in only 64 bits (`sm_key`); increasing the length beyond 32 is
possible, but would require using more than 64 bits.

Kmer length may also have an impact on the maximum number of splits. In the
*smufin* pipeline splits are subsequences of reads broken up by unknown bases
(`N`). The maximum number of splits is defined as `MAX_SPLITS` in
`src/input.hpp`, and is set to 10 by default; may need to be increased for
short kmers and inputs with an unusually high number of unknown bases.

## Length of Reads

Algorithmically speaking *smufin* doesn't take read length into consideration
and it should be able to handle reads of any length. However, certain data
structures are optimized to improve performance under different scenarios.
The `MAX_READ_LEN` environment variable, which defaults to 100, should be used
at compile time to define the maximum expected read length in the input.
