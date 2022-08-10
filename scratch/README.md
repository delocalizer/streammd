## Goal

* Construct r1+r2 fastqs from some read pairs and their dupes as marked by
  Picard MarkDuplicates to use as a test set for streammd.

## Method

* Use the `DI:i:` tag added with `TAG_DUPLICATE_SET_MEMBERS=true` to identify
  reads in the same duplicate set.

## Hurdle

* https://github.com/broadinstitute/picard/issues/1715 - MarkDuplicates
  `TAG_DUPLICATE_SET_MEMBERS` is broken for `ASSUME_SORT_ORDER=queryname`.
  So we'll have to work from a coord-sorted bam; slightly less convenient as
  the reads are no longer paired but we can work with it.

## Steps

Experiment determined that the first 2002 `DI:` tagged records of 1% sampled
coord-sorted bam had 1000 unique `DI:` tag values that occurred more than once,
suggesting 'nearby' mapped members of the duplicate set:
```bash
module load samtools/1.10
samtools view -f 1024 -s 0.01 \
	140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.coord.MD.bam|\
	grep -m 2002 -oP 'DI:i:\d+'|\
	sort -k3n,3n -t ':'|\
	uniq -c|\
	grep -v -c ' 1 '
1000
```

Create a pattern file including surrounding TAB for grepping:
```bash
samtools view -f 1024 -s 0.01\
	 140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.coord.MD.bam|\
	grep -m 2002 -oP 'DI:i:\d+'|\
	sort -k3n,3n -t':'|\
	uniq -c|\
	grep -v ' 1 '|\
	awk '{print "\t"$2"\t"}' > 1000.DI.n_gt_1
```

Experiment determined that the first 4050 matches of these DI tags formed a
complete closed set of one original pair plus at least one dupe pair (4 repeated
values of a DI: tag):
```bash
cat <(samtools view -H 140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.coord.MD.bam)\
	<(samtools view 140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.coord.MD.bam|grep -m 4050 -f 1000.DI.n_gt_1) > 1000.DI.n_gt_1.sam
```

Sort by qname again (for convenient inspection) and convert to fastq:
```bash
samtools sort -O sam -n 1000.DI.n_gt_1.sam  -T tmp -o 1000.DI.n_gt_1.qname.sam
samtools fastq -1 1000.DI.n_gt_1.1.fastq -2 1000.DI.n_gt_1.2.fastq 1000.DI.n_gt_1.qname.sam
```

Align fastqs and mark duplicates:
```bash
module load bwakit
module load picard
seqtk mergepe 1000.DI.n_gt_1.1.fastq 1000.DI.n_gt_1.2.fastq |\
	cutadapt --trim-n --interleaved \
		-a file:/reference/genomeinfo/adapters/IlluminaGenericAdapters_R1.fa \
		-A file:/reference/genomeinfo/adapters/IlluminaGenericAdapters_R2.fa -|\
	bwa mem -t 6 -p -K 100000000 \
		/reference/genomes/GRCh37_ICGC_standard_v2/indexes/BWAKIT_0.7.12/GRCh37_ICGC_standard_v2.fa - 2>/dev/null |\
	samtools sort -@ 6 -m 6G -|\
	samtools view -b - > 1000.DI.n_gt_1.coord.bam

java -jar $PICARD_HOME/picard.jar MarkDuplicates \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=1000.DI.n_gt_1.coord.bam \
	OUTPUT=1000.DI.n_gt_1.coord.MD.bam \
	METRICS_FILE=1000.DI.n_gt_1.coord.MD.bam.metrics \
	ASSUME_SORT_ORDER=coordinate \
	TAG_DUPLICATE_SET_MEMBERS=true
```
At the end of this we have:

1. r1 & r2 fastqs containing 2025 records. By construction these are all
   members of duplicate sets. When aligned, one pair will be the original
   and the other flagged 1024.
```
1000.DI.n_gt_1.1.fastq
1000.DI.n_gt_1.2.fastq
```

2. A bam aligned from these and with duplicates marked, plus dupe metrics:
```
1000.DI.n_gt_1.coord.bam
1000.DI.n_gt_1.coord.MD.bam
1000.DI.n_gt_1.coord.MD.bam.metrics
```
