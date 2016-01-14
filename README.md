# shmufin: Hash-based Mutation Finder

*shmufin* is the codename of a SMUFIN redesign based on hashtables. More
specifically, *shmufin* stands for Sparse Hash-based Mutation Finder, and is
pronounced /ʃmʌfɪn/. In turn, [SMUFIN][smufin] is a reference-free method to
identify mutations on tumor genomes, comparing them directly against the
corresponding normal genome of the same patient.

## Compile

Compiling *shmufin* requires `make`, a compiler such as `gcc` with C++11
support (>= 4.8), and the following libraries:

 - [boost::atomic][boost] (>= 1.55)
 - [sparsehash][sparsehash] (>= 2.0)
 - [folly::ConsumerProducerQueue][folly]
 - [moodycamel::ConcurrentQueue][concurrentq]

The paths for each library can be configured using a custom `make.conf` file,
see `make.conf.sample` for an example. On Debian-based systems, packages for
the former two libraries are available as:

 ```
 libboost1.55-dev
 libboost-atomic1.55.0
 libsparsehash-dev
 ```

*shmufin*'s k-mer length is determined at compile time, but can be easily
configured changing the `KMER_LEN` variable (defaults to 30).

## Run

 ```
 Usage: sm-filter [OPTIONS] -i INPUT_FILE -m MAP_FILE
 Options:
  -i, --input INPUT_FILE
  -m, --mapping MAP_FILE
  -p, --pid ID
  -l, --loaders NUM_LOADER_THREADS
  -f, --filters NUM_FILTER_THREADS
  --disable-filter
  --disable-stats
  -h, --help
 ```

## Maintainers

Jordà Polo `<jorda.polo@bsc.es>`, 2015-2016.

[smufin]: http://cg.bsc.es/smufin/ "SMUFIN"
[boost]: http://www.boost.org/ "Boost"
[sparsehash]: https://github.com/sparsehash/sparsehash "Sparse Hash"
[folly]: https://github.com/facebook/folly "Folly"
[concurrentq]: https://github.com/cameron314/concurrentqueue "ConcurrentQueue"
