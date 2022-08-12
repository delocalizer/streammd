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

* 256 (0x100) means not primary alignment (secondary, multimapper). Empirically
  we don't seem to output these from our mapping workflow.
* 1024 (0x400) is the SAM flag for 'duplicate'. MarkDuplicates TAGGING_POLICY
  optionally sets the DT attribute to LB or SQ to record as PCR or optical dupe
  respectively.
* 2048 (0x800) means supplementary (chimeric and not the representative). We do
  get lots of these in our bams.
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

## Data exploration

I extracted the first 1M dupes from a05168de-0e6d-4867-82b1-22d330dac0f8:

```
samtools view -h -f 1024 \
	/mnt/lustre/working/genomeinfo/sample/f/5/f5d4165d-0768-4e9c-8888-9810fce7210f/aligned_read_group_set/a05168de-0e6d-4867-82b1-22d330dac0f8.bam |\
	head -1000000 \
	> 1Mdups.sam
```
Now looking at the distribution of SAM flags:
```
samtools view 1Mdups.sam |cut -f2|sort|uniq -c|sort -n
```
we get:
```
      9 3137
      9 3257
     10 3201
     12 3185
     12 3209
     13 3153
     13 3169
     16 3233
     16 3249
     22 3145
     22 3155
     22 3193
     22 3219
     25 3217
     26 3171
     30 3187
     32 3251
     38 3139
     39 3235
     58 3203
    379 1201
    422 1137
    424 1089
    428 1153
    593 1209
    632 1161
   1645 1105
   1696 1185
   1750 1169
   1818 1121
   4867 1145
   5182 1097
 244290 1107
 244293 1187
 245516 1171
 245517 1123
```

So the great majority are 1123, 1171, 1187, 1107:

  * 1123 = read paired, read mapped in proper pair, mate reverse strand, first in pair, is dupe
  * 1171 = read paired, read mapped in proper pair, read reverse strand, second in pair, is dupe
  * 1187 = read paired, read mapped in proper pair, mate reverse strand, second in pair, is dupe
  * 1107 = read paired, read mapped in proper pair, read reverse strand, first in pair, is dupe

i.e. 'normal' pairs.
Then there's 1097 and 1145, which are mate unmapped.
Then there's 1121,1169,1105,1185 which are 'not proper pair' — i.e. mates mapped
too far apart or too close (lots of both cases in our data).
Then a big tail of other stuff.

# Testing our impl

20220812
```
[conradL@outgrabe streammd]$ samtools view -f 1024 scratch/1000.DI.n_gt_1.coord.MD.sam |wc -l
2034
[conradL@outgrabe streammd]$ python src/streammd/markdups.py < scratch/1000.DI.n_gt_1.coord.sam|samtools view -f 1024 -|wc -l
1931
```

So we're currently at ~ 95% of what Picard MarkDups does; 103 records short.
Let's try to find them.

Here are 4: DI:i:1327
```
samtools view -h -f16 scratch/1000.DI.n_gt_1.coord.MD.sam |samtools view -f 32|wc -l
4
```
These are 2 sets of pairs, where they're marked 0x10 (16) and 0x20 (32) — i.e.
read and its mate are both reverse. We explicitly don't count them at the
moment.

There are 0 that are both forward:
```
samtools view -h -F16 scratch/1000.DI.n_gt_1.coord.MD.sam |samtools view -F 32 -|wc -l
0
```

There are 0 that are unmapped (either end):
```
[conradL@outgrabe streammd]$ samtools view -f4 scratch/1000.DI.n_gt_1.coord.MD.sam |wc -l
0
[conradL@outgrabe streammd]$ samtools view -f8 scratch/1000.DI.n_gt_1.coord.MD.sam |wc -l
0
```

What about clipping? Are we handling that properly?
```
[conradL@outgrabe streammd]$ samtools view -f1024 scratch/1000.DI.n_gt_1.coord.MD.sam |cut -f6|grep -P '\d+[SH]'|wc -l
59
[conradL@outgrabe streammd]$ python src/streammd/markdups.py < scratch/1000.DI.n_gt_1.coord.sam |samtools view -f 1024|cut -f6|grep -P '\d+[SH]'|wc -l
21
```
So we seem a bit short there... interestingly there is also one hard clip and
two soft clips in our impl that aren't in Picard MD so that's worth looing at:
```
diff <(samtools view -f1024 scratch/1000.DI.n_gt_1.coord.MD.sam |cut -f6|grep -P '\d+[SH]'|sort|uniq) <(python src/streammd/markdups.py < scratch/1000.DI.n_gt_1.coord.sam |samtools view -f 1024|cut -f6|grep -P '\d+[SH]'|sort|uniq)|grep '^>'
> 56H45M
> 61M40S
> 98M3S
```

Look at these two pairs that are MD in Picard:
```
HWI-ST1213:151:C1DTBACXX:2:2109:18011:23563	99	chr1	59256660	60	101M	=59256870	311	ATTTATAAAGAAAAAAAAAGGTGTATTTGGCTCATAATTCTGCTGCTTGGAAAGTTCATGATTGGGCATCTGCATCTGGTAAGGGCCTCAGGCTGTTTCTA	CCCFFFFFHHHHGIJJIJIGI?FDEIJJJIJJJJJJJJJIJJJJIJIJJJEIJJGHHHHHHHFFFFFEEEECEDD5@CCCDDDDDD?BBDDCDDDDDCDD@	MD:Z:101	PG:Z:MarkDuplicates	DI:i:2271	NM:i:0	AS:i:101	DS:i:2	XS:i:31
HWI-ST1213:151:C1DTBACXX:2:2109:18011:23563	147	chr1	59256870	60	101M	=59256660	-311	TGAGAACTAATAGAGCTAGAAGAACTCACTCCCTGCCTCACCCCAGGGAGGGCATTACTCTATTCATGAAGGATCTGCTCCCATGGCCCAAACACCTCCCA	DDEEDDEDEDDCDDDDDDDDFC>5DCA<DDDDDDDFFDBBJIFJJJJJHHHIIIHHFJIJJIJJIJIIJJJJJIJJIHGIJJJJJJIIHHHHFDFFFFCCB	MD:Z:101	PG:Z:MarkDuplicates	DI:i:2271	NM:i:0	AS:i:101	DS:i:2	XS:i:31
HWI-ST1213:151:C1DTBACXX:2:2109:2472:34994	1123	chr1	59256660	60	77M24S	=59256904	311	ATTTATAAAGAAAAAAAAAGGTGTATTTGGCTCATAATTCTGCTGCTTGGAAAGTTCATGATTGGGCATCTGCATCTTGGAAAAGGGGCCCAGGCGGTTTT	B@CFFADFHGHGHJJIJJGHIBGFGIJIJJJJJJJGJJJIIIIJIIJIJJEHFD@CHHEHHEFDBDFCA@###############################	MD:Z:77	PG:Z:MarkDuplicates	DI:i:2271	NM:i:0	AS:i:77	DS:i:2	XS:i:22
HWI-ST1213:151:C1DTBACXX:2:2109:2472:34994	1171	chr1	59256904	60	34S67M	=59256660	-311	TGAAACCTAATACACCTATACACTCTCCCCCCCCGCCTCACCCCAGAGAGGGCATAACTCTATTCATGAAGCACCTGCCCCGATGCCCCAAACACCTCCCA	######################################################################################@?<2++<+D=D?B;=	MD:Z:12G8T15G1T4T2C3G15	PG:Z:MarkDuplicates	DI:i:2271NM:i:7	AS:i:32	DS:i:2	XS:i:22
```

so we kind of fixed that by using `alignment.pos - alignment.qstart + 1` but
that does not handle soft clipping on the mate. This highlights a problem with
our approach of handling the pair separately; you can end up with one read of
a pair marked as dupe and the other not:

```
[conradL@outgrabe streammd]$ samtools view -f 1024 out.sam |grep -P 'HWI-ST1213:151:C1DTBACXX:2:2109:2472:34994|HWI-ST1213:151:C1DTBACXX:2:2109:18011:23563'
HWI-ST1213:151:C1DTBACXX:2:2109:2472:34994	1171	chr1	59256904	60	34S67M	=59256660	-311	TGAAACCTAATACACCTATACACTCTCCCCCCCCGCCTCACCCCAGAGAGGGCATAACTCTATTCATGAAGCACCTGCCCCGATGCCCCAAACACCTCCCA	######################################################################################@?<2++<+D=D?B;=	NM:i:7	MD:Z:12G8T15G1T4T2C3G15	AS:i:32	XS:i:22
```

In this case only we didn't get the HWI-ST1213:151:C1DTBACXX:2:2109:2472:34994
alignment starting at chr1:59256660 marked as a duplicate even though it should
be.

So after handling groups 'properly':
```
[conradL@outgrabe streammd]$ python src/streammd/markdups.py < scratch/1000.DI.n_gt_1.qname.sam|samtools view -f 1024 -|wc -l
1978
```
now we're only 56 records short, or 97.5% of what Picard gives us.
