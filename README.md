# streammd

## Description

Single-pass duplicate marking with a Bloom filter. This is a probabilistic
approach to duplicate detection with the following features:

 * very fast (typically four to five times faster than Picard MarkDuplicates)
 * low in resources (5GB for 1B input reads)
 * tunable target false positive rate
 * streaming input and output

Detected duplicates have the SAM FLAG 0x400 bit set in the outputs.

Summary metrics are output to STDERR at the end of processing.

## Implementation

Currently only paired reads in qname-grouped order are supported. Soft
clipping is correctly handled when determining fragment ends.

## Install

```bash
git clone https://github.com/delocalizer/mumbai
cd mumbai
pip install -r requirements.txt --user .
```

## Test
(optional, requires tox)
```bash
TODO
```

## Usage

0. get help

```bash
streammd --help
```

1. mark duplicates on an input SAM file record stream 
```
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

### Handling outputs

When using the default number of consumer processes `streammd` generates
outputs to STDOUT extremely rapidly. This can mean the bottleneck in a
workflow becomes anything downstream of it — for example, if you want
to write the outputs to bam format using `samtools view` you should use
8 or more compression threads for optimal speed:
```
samtools view -h some.bam|streammd|samtools view -h -@8 -o some.MD.bam
```
