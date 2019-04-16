# Changelog

## 2.0.0-b2 -- UNRELEASED
- `count`:
  - Cached kmers are now inserted to the table as soon as the same root is
    seen again; until now, cached kmers were only inserted after finding the
    same stem, which potentially lead to missing certain kmers and forced
    stronger requirements on directional filter parameters.
  - Turn `prefilter` into a configuration option instead of a command; it can
    now be optimized and performed during table conversion.
  - Support different methods to convert tables indexed by stem to tables
    indexed by root, including: `mem`, the default in-memory conversion, as
    well as `stream` and `slice`, two new methods that dump to and read from
    disk.
- `group_rocks`:
  - Sort reads in groups, avoiding duplicates.
  - Make `dump` a explicit command instead of automatically running it.
  - Generate group sets to ease extremal group finding.

## 2.0.0-b1 -- 2018-11-23
- Use new versioning scheme and changelog.
