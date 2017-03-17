# smufin: Mutation Finder and Genome Comparison Tool

*smufin* is a mutation finder and side-by-side genome comparison and
manipulation toolset based on k-mers. It's both the current reference
implementation and a redesign of the original [SMUFIN][smufin], a
reference-free method to identify mutations on tumor genomes, comparing them
directly against the corresponding normal genome of the same patient.

*smufin* has been designed as a reconfigurable set of checkpointable stages,
and supports different modes of execution to adapt to the characteristics of
the hardware where it's running: from scale-out executions in large data
centers to scale-up solutions that take advantage of accelerators and storage
class memory in a single machine.

## Compile

Compiling *smufin* requires `make`, a compiler such as `gcc` with C++11
support (>= 4.8), and the following libraries:

 - [sparsehash][sparsehash] (>= 2.0)
 - [boost][boost] (>= 1.55): Property trees and string algorithms
 - [ConcurrentQueue][concurrentq] and [ReaderWriterQueue][rwq]: MPMC and SPSC
   queues
 - [libbf][libbf]: Bloom filters
 - [RocksDB][rocksdb] (>= 4.9): Key-value store for flash storage
 - [htslib][htslib]: Parse BAM files

The paths for each library can be configured using a custom `make.conf` file,
see `make.conf.sample` for an example. On Debian-based systems, packages for
the first two and last two libraries are available as: `libsparsehash-dev
libboost1.55-dev librocksdb-dev libhts-dev`.

*smufin*'s makefile defaults to a minimal output. For a more verbose output,
use the following:

 ```
 VERBOSE=1 make
 ```

## Run

Running *smufin* requires a configuration file such as the sample
[smufin.conf](smufin.conf). The following command line options can be used to
override the configuration file.

 ```
 Usage: sm -c CONFIG [OPTIONS]
 Options:
  -p, --partitions NUM_PARTITIONS
  --pid PARTITION_ID
  -l, --loaders NUM_LOADER_THREADS
  -s, --storers NUM_STORER_THREADS
  -f, --filters NUM_FILTER_THREADS
  -m, --mergers NUM_MERGE_THREADS
  -g, --groupers NUM_GROUP_THREADS
  --input-normal INPUT_FILES
  --input-tumor INPUT_FILES
  -o, --output OUTPUT_PATH
  -x, --exec COMMANDS
  -h, --help
 ```

Input file arguments given to `--input-{normal,tumor}` are supposed to be
either a single file, or a quoted wildcard expandable string, e.g.
`"file-*.fq.gz"` or `"file-[12].fq.gz"`.

### Commands

The argument passed to the `--exec` flag, or `core.exec` configuration option,
must be a list of stage commands separated by semicolon. Commands are
prepended by a stage name followed by colon, and chained in a comma-separated
list. E.g. `count:run,dump` or `count:restore;filter:run,dump`. The following
list contains all available stages and commands:

 * `prune`
   * `run`: generates a bloom filter of stems that have been observed in the
     input more than once; optional stage that can be run first to save memory
     during `count`.
 * `count`: build frequency table.
   * `run`: counts frequency of normal and tumoral kmers in input sequence,
     ignoring kmers whose stem is only seen once; counters hold values up to
     2^16.
   * `dump`: serialize kmer frequency as [sparsehash
     tables](doc/formats.md#sparsehash-table) indexed by stem, for
     checkpointing and/or later analysis.
   * `restore`: unserialize dumped frequency tables from disk.
   * `stats`: display frequency stats, including size of different tables, and
     histograms for normal and tumoral counts.
   * `export`: serialize frequencies as plain [CSV
     table](doc/formats.md#csv-table) files containing kmers along with normal
     and tumoral counters; rows can be limited to kmers that meet certain
     criteria through configuration options `export-{min,max}`.
 * `filter`: select breakpoint candidates and build indexes.
   * `run`: build filter normal and tumoral (mutated and non-mutated) indexes
     containing candidate reads, along with their IDs and positions of
     candidate kmers.
   * `dump`: finalize writing filter indexes to disk; when using RocksDB
     indexes, force a compaction.
   * `stats`: display sizes of the different filters.
 * `merge`: combine multiple filter indexes.
   * `run`: read and combine filter indexes from different partitions into a
     single, unified index in RocksDB. Merges all possible indexes,
     sequentially one at a time.
   * `run_{seq,k2i,i2p}_{nn,tn,tm}`: read and combine specific filter indexes
     from different partitions into a single RocksDB instance.
   * `stats`: display sizes of the merged filters.
 * `group`: match candidates that belong to the same region.
   * `run`: window-based group leader selection and retrieval of related
     reads.
   * `stats`: display number of groups generated by each thread.

Note that commands need to follow a certain order, and some stages can't be
executed without running earlier stages first. The following graph shows the
dependencies between *smufin* commands:

![Command dependency graph](doc/figures/deps.png)

## Additional Documentation

 * [Data and File Formats](doc/formats.md)
 * [Performance Tuning and Limitations](doc/performance.md)

## Maintainers

Jordà Polo `<jorda.polo@bsc.es>`, 2015-2017.

[smufin]: http://cg.bsc.es/smufin/ "SMUFIN"
[boost]: http://www.boost.org/ "Boost"
[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
[rocksdb]: https://github.com/facebook/rocksdb "RocksDB"
[concurrentq]: https://github.com/cameron314/concurrentqueue "ConcurrentQueue"
[rwq]: https://github.com/cameron314/readerwriterqueue "ReaderWriterQueue"
[libbf]: https://github.com/mavam/libbf "libbf"
[htslib]: https://github.com/samtools/htslib "htslib"
