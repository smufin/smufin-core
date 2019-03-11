# Data & File Formats

 * [Intermediate](#intermediate)
   * [Sparsehash Table](#sparsehash-table)
   * [CSV Table](#csv-table)
   * [SEQ Index](#seq-index)
   * [K2I Index](#k2i-index)
   * [I2P Index](#i2p-index)
 * [Output](#output)
   * [Groups](#groups)


## Intermediate

Intermediate files are generated during the execution of the pipeline, and can
be used for further analysis as well as checkpointing.

### Sparsehash Table

*Stage*: `count`
*Filename*: `table.<PID>-<LID>.sht`

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
   |  |  |  |  |      |  |  |  `----- [3][3][0]: T, T, Normal
   |  |  |  |  |      |  |  `-------- [3][2][1]: T, G, Tumor
   |  |  |  |  |      |  `----------- [3][2][0]: T, G, Normal
   |  |  |  |  |      `-------------- [3][1][1]: T, C, Tumor
   |  |  |  |  |
   |  |  |  |  `--------------------- [0][2][0]: A, G, Normal
   |  |  |  `------------------------ [0][1][1]: A, C, Tumor
   |  |  `--------------------------- [0][1][0]: A, C, Normal
   |  `------------------------------ [0][0][1]: A, A, Tumor
   `--------------------------------- [0][0][0]: A, A, Normal
 ```

In the following example, the kmer `ACAGGTCCAAGGAAAGTCTTAGTGTGGGGA` has a
normal counter of 2, and a tumoral counter of 1, while the kmer
`TCAGGTCCAAGGAAAGTCTTAGTGTGGGGC` has a tumoral counter of 7:

 ```
 CAGGTCCAAGGAAAGTCTTAGTGTGGGG [2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0]
 ```

### CSV Table

*Stage*: `count/export`
*Filename*: `table.<PID>-<LID>.csv`

Representation of tables as CSV files with three columns. The 1st column
represents kmers, the 2nd column normal counters, and the 3rd column tumoral
counters. E.g.:

 ```
 GCAGGTCCAAGGAAAGTCTTAGTGTGGGGT,31,1
 TGCAGGTCCAAGGAAAGTCTTAGTGTGGGG,3,29
 ```

### SEQ Index

*Stage*: `filter`, `merge`
*Filename*: `index-seq-{nn,tn,tm}.{txt.rdb}`

SEQ files map sequence IDs to sequences. E.g.

 ```
 chr20.b-20231724/2 TCTCTTCTGTGCCCTGAATTCTCTCTCTCTCCCTCTCACACACACACACACACACACACACGCACG
 ```

### K2I Index

*Stage*: `filter`, `merge`
*Filename*: `index-k2i-{nn,tn}.{txt.rdb}`

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

*Stage*: `filter`, `merge`
*Filename*: `index-i2p-tm.{txt.rdb}`

I2P stands for *ID to Positions*. I2P files contain, for each candidate read
ID, positions within the read sequence that reference candidate kmers.
For the default configuration with `MAX_READ_LEN` set to 100, positions are
encoded as four `uint64_t` bitmaps: two in direction A, and two in direction B
(hence supporting reads of up to `128+k-1` bases). E.g.:

 ```
 chr20-13310454 2047 0 1 0
 chr20-18864072 0 0 33554176 0
 ```


## Output

### Groups

*Stage*: `group`
*Filename*: `group.<PID>-<GID>.json`

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

*Stage*: `group_rocks`
*Filename*: `group.<PID>-<GID>.msgpack`

A MSGPACK file containing a list of groups. The schema of each MSGPACK group
is as follows:

```
"lead" : READ,
"pos" : [ POS-A, POS-B ],
"kmers" : [ KMERS-A, KMERS-B ],
"reads" : [ READS-N, READS-T ]
```

Where upper-case words stand for:

 - READ: Tuple of two strings representing a read ID and its sequence,
 - POS-[A|B]: Lists of integers corresponding to candidate positions in a
   sequence, in direction A and B,
 - KMERS-[A|B]: Lists of KMER in direction A and B,
 - READS-[N|T]: Lists of READ, normal and tumoral.

*Stage*: `group_rocks`
*Filename*: `sets.<PID>-<GID>.txt`

A sets file is a text file containing a summary of the information available
in group files, and its meant to ease discarding groups that aren't relevant,
e.g. groups whose reads are all contained in another group. Each line of the
file represents a group, and it has the following format:

```
<lead ID> <number of items> <space-separated list read IDs>
```

[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
