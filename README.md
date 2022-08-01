## Goal

Streaming duplicate marking with a Bloom filter. This is a probabilistic
approach to duplicate detection but it should be very fast and lower mem. Also,
the target FP rate is tunable to whatever you wish so we should be able to make
something "good enough".

## Back of the envelope

60x human WGS 2x150bp reads ~ 6e8 fragments, 1.2e9 se reads

## Notes

* 1024 is the SAM flag for 'duplicate'. Picard MarkDuplicates TAGGING_POLICY
  optionally sets the DT attribute to LB or SQ to record as PCR or optical dupe
  respectively.

An example:

    samtools view \
        /mnt/lustre/working/genomeinfo/sample/f/5/f5d4165d-0768-4e9c-8888-9810fce7210f/aligned_read_group_set/a05168de-0e6d-4867-82b1-22d330dac0f8.bam|\
        grep -m 4 -P 'D8QSB6V1:72:C1K64ACXX:4:1211:8635:39730|D8QSB6V1:72:C1K64ACXX:4:2206:19121:40710'

## Implementation thoughts

* multiple reader threads with the BF in shared mem

## Questions

* Require sorted input/output?

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

