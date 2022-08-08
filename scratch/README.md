# Goal

* Construct r1+r2 fastqs from some read pairs and their dupes as marked by
  Picard MarkDuplicates to use as a test set for streammd.

# Method

* Use the `DI:i:` tag added with `TAG_DUPLICATE_SET_MEMBERS=true` to identify
  reads in the same duplicate set.

# Hurdle

* https://github.com/broadinstitute/picard/issues/1715 - MarkDuplicates
  `TAG_DUPLICATE_SET_MEMBERS` is broken for `ASSUME_SORT_ORDER=queryname`.
