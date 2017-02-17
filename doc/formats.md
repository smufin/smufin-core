# Data & File Formats

## Input

### Input File

Input files consist of a newline-seperated list files in one of the supported
formats: gzipped FASTQ files or indexed BAM files. There are two different
kinds of input files: normal and tumoral. Normal files contain the string
`_N_` in their name, while tumoral files contain `_T_`. E.g. sample input file
with FASTQ samples:

 ```
 ./test/00_N_insertion.fq.gz
 ./test/00_T_insertion.fq.gz
 ```

## Intermediate

### Sparsehash Table

*Stage*: `count`
*Filename*: `table.PID-LID.sht`

Tables contain normal and tumoral kmer frequencies, and can be dumped/restored
using the default [sparsehash][sparsehash] serialization.

Kmers are indexed using their stem as key, so each entry in the table includes
values for all possible kmers that can be derived from a stem, for a total of
32 counters. Counters are indexed as 3-dimensional `4x4x2` arrays, where the
first and second dimension stand for the numerical code of the initial and
last base in the kmer, respectively, and the third dimension represents the
two kinds of counters, normal and tumoral. Bases are encoded as `A: 0, C: 1,
G: 2, T: 3`, while kinds are encoded as `N: 0, T: 1`, as follows:

 ```
 [00,01,02,03,04,...,27,28,29,30,31]
   |  |  |  |  |      |  |  |  |  |
   |  |  |  |  |      |  |  |  |  `-- [3][3][1]: T, T, Tumor
   |  |  |  |  |      |  |  |  `----- [3][2][0]: T, G, Normal
   |  |  |  |  |      |  |  `-------- [3][1][1]: T, C, Tumor
   |  |  |  |  |      |  `----------- [3][0][0]: T, A, Normal
   |  |  |  |  |      `-------------- [2][3][1]: G, T, Tumor
   |  |  |  |  |
   |  |  |  |  `--------------------- [0][1][0]: C, A, Normal
   |  |  |  `------------------------ [0][3][1]: A, T, Tumor
   |  |  `--------------------------- [0][2][0]: A, G, Normal
   |  `------------------------------ [0][1][1]: A, C, Tumor
   `--------------------------------- [0][0][0]: A, A, Normal
 ```

In the following example, the kmer `ACAGGTCCAAGGAAAGTCTTAGTGTGGGGA` has a
normal counter of 2, and a tumoral counter of 1:

 ```
 CAGGTCCAAGGAAAGTCTTAGTGTGGGG [2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
 ```

### CSV Table

*Stages*: `count/export`
*Filenames*: `table.PID-LID.csv`

Stores tables as CSV files with three columns. The 1st column represents
kmers, the 2nd column normal counters, and the 3rd column tumoral counters.
E.g.:

 ```
 GCAGGTCCAAGGAAAGTCTTAGTGTGGGGT,31,1
 TGCAGGTCCAAGGAAAGTCTTAGTGTGGGG,3,29
 ```

### K2I Index

*Stages*: `filter`, `merge`
*Filenames*: `filter-k2i-{nn,tn}.{txt.rdb}`

K2I stands for *Kmer to IDs*. K2I files contain, for each candidate kmer, a
list of sequence IDs that contain that kmer. Every line in a K2I file
represents a candidate kmer and must have at least 2 columns: the first column
is the kmer itself, and the second column is the number of associated reads
that contain the kmer. The remaining columns are a variable space-separated
list of sequence IDs. E.g.:

 ```
 GGGGTGCAGGTCCAAGGAAAGTCTTAGTGT 3 chr20-18462176 chr20-5273350 chr20-6534694
 TGGGGGTGCAGGTCCAAGGAAAGTCTTAGT 1 chr20-18462176
 ```

### I2P Index

*Stages*: `filter`, `merge`
*Filenames*: `filter-i2p-tm.{txt.rdb}`

I2P stands for *ID to Positions*. I2P files contain, for each candidate read
ID, positions within the read sequence that reference candidate kmers.
Positions are encoded as four `uint64_t` bitmaps: two in direction A, and two
in direction B (hence supporting reads of up to `128+k-1` bases). E.g.:

 ```
 chr20-13310454 2047 0 1 0
 chr20-18864072 0 0 33554176 0
 ```


## Output

### Groups File

*Stages*: `group`
*Filenames*: `group.PID-GID.json`

A JSON file containing a dict of groups indexed with the sequence ID of their
leader as key. The schema of a JSON groups file is as follows:

```
{
    "ID" : {
        "lead" : READ,
        "pos-a" : POS,
        "pos-b" : POS,
        "kmers-a" : [ KMER, KMER, ... ],
        "kmers-b" : [ KMER, KMER, ... ],
        "reads-n" : [ READ, READ, ... ],
        "reads-t" : [ READ, READ, ... ]
    },
    "ID" : { ... },
    ...
}
```

Where upper-case words stand for:

 - ID: Unique string that identifies a sequence, e.g. `chr20.b-20231724/2`,
 - POS: List of integers corresponding to candidate positions in a sequence,
   e.g. `[ 19, 20, 21 ]`.
 - KMER: Tuple formed by a string and four integers. The string of length K
   represents a candidate kmer. The first two integers represent number of
   retrieved reads (normal and tumoral), while the other two represent
   discarded reads. E.g. `[ "CTCTCCCTCTCACACACACACACACACA", 20, 19, 0, 0 ]`.
 - READ: Tuple of 2 strings: ID and read sequence; e.g.
```
[
   "chr20.b-20231724/2",
   "TCTCTTCTGTGCCCTGAATTCTCTCTCTCTCCCTCTCACACACACACACACACACACACACGCACG",
]
```

[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
