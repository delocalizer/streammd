# streammd

Single-pass probabilistic duplicate marking of alignments with a Bloom filter.

Input is a SAM file format input stream: a valid header followed by reads in
qname-grouped order (e.g. output of `bwa mem`). Detected duplicates have the
SAM FLAG 0x400 bit set in the outputs and summary metrics are written to STDERR
at the end of processing.

### Features

* Very fast (when run with default settings `streammd` is typically 3-4x faster
  than Picard MarkDuplicates).
* Low memory use, especially for large libraries (approx 5GB for 1B read pairs).
* Soft clipping is correctly handled when determining ends.
* Tunable target false positive rate.
* Streaming input and output.

### Limitations

Inherent, due to the nature of the single-pass operation:

* `streammd` retains the first encountered template as the original and marks
  subsequently encountered copies as duplicates. This differs from Picard
  MarkDuplicates which marks those of lowest quality as the duplicates.

Implementation specific:

* Output is not deterministic when using more than 1 consumer process.
* Only paired reads are supported.

## Install

```bash
git clone https://github.com/delocalizer/mumbai
cd mumbai
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

|    n    |    p    |    mem    |
| ------- | ------- | --------- |
|1.00E+06 |1.00E-06 |0.0033GB   |
|1.00E+06 |1.00E-07 |0.0039GB   |
|1.00E+06 |1.00E-08 |0.0045GB   |
|1.00E+06 |1.00E-09 |0.0050GB   |
|1.00E+07 |1.00E-06 |0.0335GB   |
|1.00E+07 |1.00E-07 |0.0391GB   |
|1.00E+07 |1.00E-08 |0.0446GB   |
|1.00E+07 |1.00E-09 |0.0502GB   |
|1.00E+08 |1.00E-06 |0.3348GB   |
|1.00E+08 |1.00E-07 |0.3905GB   |
|1.00E+08 |1.00E-08 |0.4463GB   |
|1.00E+08 |1.00E-09 |0.5021GB   |
|1.00E+09 |1.00E-06 |3.3475GB   |
|1.00E+09 |1.00E-07 |3.9055GB   |
|1.00E+09 |1.00E-08 |4.4634GB   |
|1.00E+09 |1.00E-09 |5.0213GB   |

As a guide, 60x human WGS 2x150bp paired-end sequencing consists of
n &#8776; 6.00E+08 templates.

### Handling outputs

When using the default settings `streammd` writes outputs to STDOUT extremely
rapidly. In a pipeline this can mean anything downstream will become a
bottleneck if it can't process its inputs fast enough — for example, if you
want to write the outputs to bam format using `samtools view` you should use 8
or more compression threads for optimal speed:

```bash
samtools view -h some.bam|streammd|samtools view -h -@8 -o some.MD.bam
```
