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
   * `dump`: write filter indexes to disk
   * `stats`: display sizes of the different filters.
 * `merge`
   * `run`: read and combine filters from different partitions into a single,
     unified filter in a RocksDB database.
 * `group`
   * `run`: window-based group leader selection and retrieval of related
     reads.

## Maintainers

Jordà Polo `<jorda.polo@bsc.es>`, 2015-2016.

[smufin]: http://cg.bsc.es/smufin/ "SMUFIN"
[boost]: http://www.boost.org/ "Boost"
[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
[folly]: https://github.com/facebook/folly "Folly"
[rocksdb]: https://github.com/facebook/rocksdb "RocksDB"
[concurrentq]: https://github.com/cameron314/concurrentqueue "ConcurrentQueue"
