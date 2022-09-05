# streammd

Single-pass probabilistic duplicate marking of alignments with a Bloom filter.

Input is a SAM file format input stream: a valid header followed by reads in
qname-grouped order (e.g. output of `bwa mem`). Detected duplicates have the
SAM FLAG 0x400 bit set in the outputs and summary metrics are written to a file
at the end of processing.

### Features

* Very fast — with default settings `streammd` is more than 2x as fast as
  Picard MarkDuplicates, and much faster is easily achievable using more
  workers.
* Low memory use, even for large libraries — with default settings `streammd`
  requires less than 4GB to process 1B templates.
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

* Output is not deterministic when using more than 1 worker process.

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
samtools view -h some.bam|streammd --paired
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
n &#8776; 6.00E+08 templates. Run the included `memcalc` tool to get an
estimate of `streammd` memory use for a given `n` and `p`.

### Pipelining

By using multiple worker processes `streammd` is capable of generating outputs
at a very high rate. For efficient pipelining, downstream tools should be run
with sufficient cpu resources to handle their inputs — for example if you want
to write the outputs to bam format using `samtools view` you should specify
extra compression threads for optimal throughput:

```bash
samtools view -h some.bam|streammd --paired|samtools view -@2 -o some.MD.bam
```

By the same token, there's no value in increasing `streammd` throughput by
specifying larger values for `-w, --workers` if a downstream tool does not have
the processing speed to handle it.
