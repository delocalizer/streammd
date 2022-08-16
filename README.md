# streammd

## Description

Single-pass duplicate marking with a Bloom filter. This is a probabilistic
approach to duplicate detection with the following features:

 * very fast (typically 2-3x faster than Picard MarkDuplicates)
 * low memory use (approx 5GB for 300M read pairs)
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

For calculating `n`, multiply the number of read pairs in the input by 3
because for each template 3 items are stored: 2 x read ends and 1 x template
ends.

For example; 60x human WGS 2x150bp paired-end sequencing ~> 6e8 templates ~> n = 1.8e9

### Handling outputs

When using the default number of consumer processes `streammd` writes outputs
to STDOUT extremely rapidly. In a pipeline this can mean anything downstream
will become a bottleneck if it can't consume the inputs fast enough — for
example, if you want to write the outputs to bam format using `samtools view`
you should use 8 or more compression threads for optimal speed:

```bash
samtools view -h some.bam|streammd|samtools view -h -@8 -o some.MD.bam
```
