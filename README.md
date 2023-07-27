# streammd

Single-pass probabilistic duplicate marking of alignments with a Bloom filter.

Input is a SAM file format input stream: a valid header followed by reads in
qname-grouped order (e.g. output of `bwa mem`). Detected duplicates have the
SAM FLAG 0x400 bit set and PG:Z: tag added, and summary metrics are written
to a file at the end of processing.

### Features

* Fast — with default settings `streammd` is ~5x faster than Picard
  MarkDuplicates.
* Low memory footprint even for large libraries — with default settings
  `streammd` requires just 4G to process 1B templates.
* High concordance with Picard MarkDuplicates metrics.
* Soft-clipped reads are correctly handled.
* Tunable target false positive rate.
* Streaming input and output.

### Limitations

Due to the nature of the single-pass operation:

* `streammd` retains the first encountered template as the original and marks
  subsequently encountered copies as duplicates. This differs from Picard
  MarkDuplicates which marks those of lowest quality as the duplicates.
* `streammd` does not differentiate between optical duplicates and PCR
  duplicates.

Additionally, the current implementation handles SAM format input only.

## Requirements

### Development

* autotools

### Installation

* Requires c++17  with `<charconv>`. gcc >= 8.1 or clang >= 7 should work.

## Install

Download a distribution tarball `streammd-[x.y.z].tar.gz` from
https://github.com/delocalizer/streammd/releases


```bash
tar -xzvf streammd-[x.y.z].tar.gz
cd streammd-[x.y.z]
./configure
# build
make
# run unit tests (optional)
make check
# install
sudo make install
```

The configure script and makefiles are generated by autotools and should
be quite portable. The usual variables and options may be used to customize
the build and install steps e.g.

```
./configure --prefix=$HOME/.local
make CXX='clang++' CXXFLAGS='-g -O3'
```

## Usage

0. get help

```bash
streammd --help
```

```
Usage: streammd [-h] [--input INPUT] [--output OUTPUT] [--fp-rate FP_RATE] [--mem MEM]
                [--allow-overcapacity] [--metrics METRICS_FILE] [--remove-duplicates]
                [--show-capacity] [--single] [--strip-previous]

Read a SAM file from STDIN, mark duplicates in a single pass and stream processed records to STDOUT.
Input must begin with a valid SAM header followed by qname-grouped records. Default log level is
'info' — set to something else (e.g. 'debug') via SPDLOG_LEVEL environment variable.

Optional arguments:
  -h, --help            	shows help message and exits 
  -v, --version         	prints version information and exits 
  --input INPUT         	Input file. [default: STDIN] 
  --output OUTPUT       	Output file. [default: STDOUT] 
  -p, --fp-rate FP_RATE 	The maximum acceptable marginal false-positive rate. [default: 1e-06]
  -m, --mem MEM         	Memory allowance for the Bloom filter, e.g "4GiB". Both binary
                                (kiB|MiB|GiB) and decimal (kB|MB|GB) formats are understood. As a
                                result of implementation details, a value that is an exact power
                                of 2 (512MiB, 1GiB, 2GiB etc) gives a modest processing speed
                                advantage (~5%) over neighbouring values. [default: "4GiB"]
  --allow-overcapacity  	Warn instead of error when Bloom filter capacity is exceeded.
                                [default: false] 
  --metrics METRICS_FILE	Output metrics file. [default: "streammd-metrics.json"]
  --remove-duplicates   	Omit detected duplicates from the output. 
  --show-capacity       	Do no work, just print the capacity of the Bloom filter that would
                                be constructed with the given --fp-rate and --mem values. 
  --single              	Accept single-ended reads as input. [default: paired-end] 
  --strip-previous      	Unset duplicate flag for any reads that have it set and are no
                                longer considered duplicate. Only ever required if records have
                                previously been through a duplicate marking step. [default: false]
```

1. mark duplicates on an input SAM record stream 

```bash
bwa mem ref.fa r1.fq r2.fq|streammd
```

## Notes

### Bloom filter capacity.

Bloom filter capacity `n` (in this context, the number of templates that can be
effectively processed) depends on the desired maximum false positive rate `p`
as well as the memory allowance. For a given `p` you should set the memory to
a value sufficient to handle the expected number of templates. Run the tool
with `--show-capacity` to view Bloom filter capacity for any `-p` and `-m`
values.

|    p     |   mem  |    n     |
| -------- | ------ | -------- |
| 1.00E-02 | 128MiB | 1.07E+08 |
| 1.00E-04 | 256MiB | 1.09E+08 |
| 1.00E-06 | 512MiB | 1.24E+08 |
| 1.00E-02 | 1GiB   | 8.56E+08 |
| 1.00E-04 | 2GiB   | 8.72E+08 |
| 1.00E-06 | 4GiB   | 9.94E+08 |

For example, as a guide: 60x human WGS 2x150bp paired-end sequencing consists
of n &#8776; 6.00E+08 templates, and the default memory setting of 4GiB is
sufficient to process this at the default false positive rate of 1.00E-06.

### False-positive rate.

`-p, --fp-rate` sets the value for the *maximum acceptable marginal
false-positive rate*. This is the value at which the tool raises an error when
with `n` items already stored, the predicted FP probability of a newly tested
item exceeds that value. This is not the same as the true bulk FP rate in all
items processed by the Bloom filter, which will always be lower as long as we
stop adding items when the marginal rate is exceeded.

The `--allow-overcapacity` option overrides this default error behaviour,
generating only a warning when the target maximum marginal FP rate is exceeded.
This is not recommended for general use.

### Hash functions.

The current implementation of `streammd` uses a fixed number of hash functions
(`k=10`) regardless of the value of `p` (`--fp-rate`) or `m` (`--mem`). Other
implementations are certainly possible, such as allowing a free choice of `k`,
or using the value of `k` that maximizes capacity `n` for a given `m` and `p`
(i.e. the capacity-optimal value). Each of these has its own drawbacks. In the
case of a free choice, the user must understand and experiment with the
non-trivial relationship between `k` and filter capacity, false-positive rate,
memory use, and performance to pick a sensible value. Most users probably just 
want reasonable performance and minimal parameter space searching, given their
computational constraint (mem) and analytic goal (FPR).

Moreover, the value of `k` that maximizes `n` for a given `p` and `m` is not
a great default for performance reasons because Bloom filter capacity is only
weakly dependent on `k` close to the capacity-optimal value, but computational
cost is linear in `k`. For example, with defaults of `m=4GiB` and `p=1e-6`,
the capacity-optimal Bloom filter requires `k=20` to store `1.2e9` items, but
halving the number of hashes to `k=10` (doubling the speed) reduces capacity by
just 17% to `1e9` items.

In summary, a fixed `k=10` is a pragmatic compromise between performance and
memory given values of `p` and `n` likely to be useful for large short-read
datasets.

### Metrics

* `TEMPLATES` = total count of templates seen.
* `TEMPLATES_UNMAPPED` = count of templates containing no mapped read.
* `TEMPLATES_MARKED_DUPLICATE` = count of templates marked duplicate.
* `TEMPLATE_DUPLICATE_FRACTION` = `TEMPLATES_MARKED_DUPLICATE / (TEMPLATES - TEMPLATES_UNMAPPED)`

### Pipelining

`streammd` is capable of generating outputs at a very high rate. For efficient
pipelining, downstream tools should be run with sufficient cpu resources to
handle their inputs — for example if you want to write the outputs to BAM
format using `samtools view` you should specify extra compression threads for
optimal throughput; for example:

```bash
bwa mem ref.fa r1.fq r2.fq|streammd|samtools view -@2 -o some.MD.bam
```
### Citing

Please consider citing https://doi.org/10.1093/bioinformatics/btad181 when referring to `streammd` in publications.
