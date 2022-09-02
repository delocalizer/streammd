# streammd

Single-pass probabilistic duplicate marking of alignments with a Bloom filter.

Input is a SAM file format input stream: a valid header followed by reads in
qname-grouped order (e.g. output of `bwa mem`). Detected duplicates have the
SAM FLAG 0x400 bit set in the outputs and summary metrics are written to a file
at the end of processing.

### Features

* Very fast — with default settings (2 hashing processes) `streammd` is
  typically &#8776; 3x faster than Picard MarkDuplicates; even faster is easily
  achievable using more hashers.
* Low memory use, even for large libraries — with default settings (_p_=10<sup>-6</sup>)
  `streammd` requires less than 4GB to process 1B templates.
* High concordance with Picard MarkDuplicates metrics.
* Soft-clipped reads are correctly handled.
* Tunable target false positive rate.
* Streaming input and output.

### Limitations

Inherent, due to the nature of the single-pass operation:

* `streammd` retains the first encountered template as the original and marks
  subsequently encountered copies as duplicates. This differs from Picard
  MarkDuplicates which marks those of lowest quality as the duplicates.
* `streammd` does not differentiate between optical duplicates and PCR
  duplicates.

Implementation specific:

* Output is not deterministic when using more than 1 consumer process.

## Install

```bash
git clone https://github.com/delocalizer/streammd
cd streammd
pip install -r requirements.txt --user .
```

## Test

(optional, requires tox)
```bash
tox
```

## Usage

0. get help

```bash
streammd --help
```

1. mark duplicates on an input SAM file record stream 

```bash
samtools view -h some.bam|streammd
```

## Notes

### Memory usage

Memory usage depends on the number of items to be stored `n` and target
maximum false positive rate `p`:

|    n     |    p     |   mem    |
| -------- | -------- | -------- |
| 1.00E+07 | 1.00E-03 | 0.017GB  |
| 1.00E+07 | 1.00E-06 | 0.033GB  |
| 1.00E+07 | 1.00E-09 | 0.050GB  |
| 1.00E+08 | 1.00E-03 | 0.167GB  |
| 1.00E+08 | 1.00E-06 | 0.335GB  |
| 1.00E+08 | 1.00E-09 | 0.502GB  |
| 1.00E+09 | 1.00E-03 | 1.674GB  |
| 1.00E+09 | 1.00E-06 | 3.348GB  |
| 1.00E+09 | 1.00E-09 | 5.021GB  |


As a guide, 60x human WGS 2x150bp paired-end sequencing consists of
n &#8776; 6.00E+08 templates.

### Handling outputs

When using the default settings `streammd` writes outputs to STDOUT extremely
rapidly. In a pipeline this can mean anything downstream will become a
bottleneck if it can't process its inputs fast enough — for example, if you
want to write the outputs to bam format using `samtools view` you should use 
extra compression threads for optimal throughput:

```bash
samtools view -h some.bam|streammd|samtools view -@4 -o some.MD.bam
```
