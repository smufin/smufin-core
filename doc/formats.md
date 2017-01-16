# File Formats

## Input Files

Input files consist of a newline-seperated list of gzipped FASTQ files. There
are two different kinds of FASTQ files: normal and tumoral. Normal files
contain the string `_N_` in their name, while tumoral files contain `_T_`.
E.g. sample input file:

 ```
 ./test/00_N_insertion.fq.gz
 ./test/00_T_insertion.fq.gz
 ```

## I2P Files

I2P stands for *ID to Positions*. I2P files contain, for each candidate read
ID, positions within the read sequence that reference candidate kmers.
Positions are encoded as four `uint64_t` bitmaps: two in direction A, and two
in direction B (hence supporting reads of up to `128+k-1` bases). E.g.:

 ```
 chr20-13310454 2047 0 1 0
 chr20-18864072 0 0 33554176 0
 ```

## K2I Files

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

## Groups Files

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
