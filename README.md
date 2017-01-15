# shmufin: Hash-based Mutation Finder

*shmufin* is the codename of a SMUFIN redesign based on hashtables. More
specifically, *shmufin* stands for Sparse Hash-based Mutation Finder, and is
pronounced /ʃmʌfɪn/. In turn, [SMUFIN][smufin] is a reference-free method to
identify mutations on tumor genomes, comparing them directly against the
corresponding normal genome of the same patient.

## Compile

Compiling *shmufin* requires `make`, a compiler such as `gcc` with C++11
support (>= 4.8), and the following libraries:

 - [sparsehash][sparsehash] (>= 2.0)
 - [boost][boost] (>= 1.55): Property trees and string algorithms
 - [moodycamel::ConcurrentQueue][concurrentq]: MPMC queue
 - [folly::ProducerConsumerQueue][folly]: SPSC queue
 - [RocksDB][rocksdb] (>= 4.9): Key-value store for flash storage

The paths for each library can be configured using a custom `make.conf` file,
see `make.conf.sample` for an example. On Debian-based systems, packages for
the former two and last libraries are available as:

 ```
 libsparsehash-dev
 libboost1.55-dev
 librocksdb-dev
 ```

*shmufin*'s makefile defaults to a minimal output. For a more verbose output,
use the following:

 ```
 VERBOSE=1 make
 ```

## Run

 ```
 Usage: sm -c CONFIG [OPTIONS]
 Options:
  -p, --partitions NUM_PARTITIONS
  --pid PARTITION_ID
  -l, --loaders NUM_LOADER_THREADS
  -s, --storers NUM_STORER_THREADS
  -f, --filters NUM_FILTER_THREADS
  -m, --mergers NUM_MERGE_THREADS
  -i, --input INPUT_FILE
  -o, --output OUTPUT_PATH
  -x, --exec COMMANDS
  -h, --help
 ```

### Commands

 * `count`
   * `run`: counts frequency of normal and tumoral kmers in input sequence,
     ignoring kmers whose stem is only seen once; counters hold values up to
     2^16.
   * `dump`: serialize hashtables that contain kmer frequencies to the
     `core.output` directory; filenames have the following format:
     `table.<PID>-<SID>.sht`, where `PID` stands for partition ID, and `SID`
     for storer ID.
   * `restore`: unserialize hashtables from disk.
   * `stats`: display frequency stats, including size of different tables, and
     histograms for normal and tumoral counts.
 * `filter`
   * `run`: build filter normal and tumoral (mutated and non-mutated) indexes
     containing candidate reads, along with their IDs and positions of
     candidate kmers.
   * `dump`: finalize writing filter indexes to disk; when using RocksDB
     indexes, force a compaction.
   * `stats`: display sizes of the different filters.
 * `merge`
   * `run`: read and combine filter indexes from different partitions into a
     single, unified index in RocksDB. Merges all possible indexes,
     sequentially one at a time.
   * `run_{seq,k2i,i2p}_{nn,tn,tm}`: read and combine specific filter indexes
     from different partitions into a single RocksDB instance.
   * `stats`: display sizes of the merged filters.
 * `group`
   * `run`: window-based group leader selection and retrieval of related
     reads.

## File Formats

### Input Files

Input files consist of a newline-seperated list of gzipped FASTQ files. There
are two different kinds of FASTQ files: normal and tumoral. Normal files
contain the string `_N_` in their name, while tumoral files contain `_T_`.
E.g. sample input file:

 ```
 ./test/00_N_insertion.fq.gz
 ./test/00_T_insertion.fq.gz
 ```

### I2P Files

I2P stands for *ID to Positions*. I2P files contain, for each candidate read
ID, positions within the read sequence that reference candidate kmers.
Positions are encoded as four `uint64_t` bitmaps: two in direction A, and two
in direction B (hence supporting reads of up to `128+k-1` bases). E.g.:

 ```
 chr20-13310454 2047 0 1 0
 chr20-18864072 0 0 33554176 0
 ```

### K2I Files

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

### Groups File

A JSON file containing a dict of groups indexed with the sequence ID of their
leader as key. The schema of a JSON groups file is as follows:

```
{
    "ID" : {
        "lead" : READ,
        "pos-f" : POS,
        "pos-r" : POS,
        "kmers-f" : [ KMER, KMER, ... ],
        "kmers-r" : [ KMER, KMER, ... ],
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

## Maintainers

Jordà Polo `<jorda.polo@bsc.es>`, 2015-2017.

[smufin]: http://cg.bsc.es/smufin/ "SMUFIN"
[boost]: http://www.boost.org/ "Boost"
[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
[folly]: https://github.com/facebook/folly "Folly"
[rocksdb]: https://github.com/facebook/rocksdb "RocksDB"
[concurrentq]: https://github.com/cameron314/concurrentqueue "ConcurrentQueue"
