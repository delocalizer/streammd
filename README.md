## Goal

Streaming duplicate marking with a Bloom filter. This is a probabilistic
approach to duplicate detection but it should be very fast and lower mem. Also,
the target FP rate is tunable to whatever you wish so we should be able to make
something "good enough".

## Back of the envelope

60x human WGS 2x150bp reads ~ 6e8 fragments, 1.2e9 se reads

## Implementation thoughts

* multiple reader threads with the BF in shared mem

## Questions

* Require sorted input/output?

* Do we need to distinguish optical dupes (DT:SQ) from PCR (DT:LB) ?
  i.e. do we use library complexity? Also, how many optical dupes do we get on
  modern SR sequencing platforms anyway? Need to do some experiments to check.

* Does Picard's implementation allow mismatches in sequence? If yes that seems
  like it would be hard to replicate with a Bloom filter.

* More generally, what exactly is the Picard MD algorithm? How does it work for
  pe reads?
