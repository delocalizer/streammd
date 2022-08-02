## Goal

Streaming duplicate marking with a Bloom filter. This is a probabilistic
approach to duplicate detection but it should be very fast and lower mem. Also,
the target FP rate is tunable to whatever you wish so we should be able to make
something "good enough".

## Back of the envelope

60x human WGS 2x150bp pe sequencing ~> 6e8 fragments (templates), 1.2e9 reads

Probability of false positive for Bloom filter with n elems already inserted:

    n: # elements inserted
    m: # bits in the array
    k: # hash functions
    pr(FP) ~ (1 - e**(-kn/m))**k

One can derive optimal values for m/n (bits per element) and k (# hash funcs):

    E: pr(FP) i.e. desired false positive probability
    k_optimal     = -log_2(E)           # optimal number of hash functions
    (m/n)_optimal ~ -1.44 log_2(E)      # optimal bits per element

If we take:

    n = 1e9  (fragments)
    E = 1e-6 (desired error rate)
    => k = 20, m = 3.3GB

    n = 1e9  (fragments)
    E = 1e-7 (desired error rate)
    => k = 24, m = 3.9GB

    n = 1e9  (fragments)
    E = 1e-8 (desired error rate)
    => k = 27, m = 4.5GB

    n = 1e9  (fragments)
    E = 1e-9 (desired error rate)
    => k = 30, m = 5GB

## Notes

* 1024 is the SAM flag for 'duplicate'. Picard MarkDuplicates TAGGING_POLICY
  optionally sets the DT attribute to LB or SQ to record as PCR or optical dupe
  respectively.
* bam @HD VN: SO: is added by samtools sort
* a05168de-0e6d-4867-82b1-22d330dac0f8.bam is from 2x100bp pe seq of Melanoma
  control cell line.
  It has:
  -        120GB size
  -    825161759 unique QNAME (templates)
  -   1654725103 alignments
  -   1506118083 alignments not marked duplicate (0x400 bit unset)
  -    148607020 alignments marked as duplicate (0x400 bit set)
  -      4226341 QNAME with > 2 alignments (supplementary + secondary)
  - 146080696649 bases (from ngscheck)

An example:

    samtools view \
        /mnt/lustre/working/genomeinfo/sample/f/5/f5d4165d-0768-4e9c-8888-9810fce7210f/aligned_read_group_set/a05168de-0e6d-4867-82b1-22d330dac0f8.bam|\
        grep -m 4 -P 'D8QSB6V1:72:C1K64ACXX:4:1211:8635:39730|D8QSB6V1:72:C1K64ACXX:4:2206:19121:40710'

## Implementation thoughts

* multiple reader threads with the BF in shared mem
* https://github.com/sangupta/bloomfilter

## Questions

* Require sorted input/output? QNAME sort will be far easier as you don't have
  to wait until the next read in a pair arrives (might be on another chrom!).
  This means we should use the unsorted output from bwa mem, which is the same
  order as the input fastqs i.e. query grouped.

* Do we need to distinguish optical dupes (DT:SQ) from PCR (DT:LB) ?
  i.e. do we use library complexity? Also, how many optical dupes do we get on
  modern SR sequencing platforms anyway? Need to do some experiments to check.

* Does Picard's implementation allow mismatches in sequence? If yes that seems
  like it would be hard to replicate with a BF.

* Picard's default DUPLICATE_SCORING_STRATEGY is SUM_OF_BASE_QUALITIES (highest
  qual is the 'original', others are the dupes) but most obvious BF impl would
  be RANDOM. How much of a difference does it make, and given that, do we care?

* More generally, what exactly is the Picard MD algorithm? How does it work for
  pe reads? From MarkDuplicates.java src:
  > The MarkDuplicates tool works by comparing sequences in the 5 prime
  > positions of both reads and read-pairs in a SAM/BAM file. [...]
  > Note that this is different from directly checking if the sequences match,
  > which MarkDuplicates does not do.
  That's a bit confusing. It mentions comparing sequence, then seems to say
  it doesn't compare sequence. ReadEndsForMarkDuplicates structure certainly
  doesn't store any sequence info.
  read1Coordinate = POS
  read1ReferenceIndex = RNAME
  clarification from http://broadinstitute.github.io/picard/faq.html:
  > Q: How does MarkDuplicates work?
  > A: The MarkDuplicates tool finds the 5' coordinates and mapping
  > orientations of each read pair (single-end data is also handled). It takes
  > all clipping into account as well as any gaps or jumps in the alignment.
  > Matches all read pairs using the 5' coordinates and their orientations.
  > It marks all but the "best" pair as duplicates, where "best" is defined as
  > the read pair having the highest sum of base qualities of bases with Q ≥ 15.
  so it just goes off end coords and orientations

