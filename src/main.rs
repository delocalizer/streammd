use std::fs::{File, OpenOptions};
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use bytesize::ByteSize;
use clap::Parser;
use const_format::formatcp;
use env_logger::Env;
use log::{error, info, warn};
use serde_json::{json, to_string_pretty};

use bloomfilter::BloomFilter;

mod markdups;

type Reader = Box<dyn BufRead>;
type Writer = Box<dyn Write>;

const DEFAULT_K: u8 = 10;
const PKG_NAME: &str = env!("CARGO_PKG_NAME");

/// Read a SAM file from STDIN, mark duplicates in a single pass and stream
/// processed records to STDOUT. Input must begin with a valid SAM header
/// followed by qname-grouped records. Default log level is 'info' — set to
/// something else (e.g. 'debug') via the RUST_LOG environment variable.
#[derive(Parser, Debug)]
#[command(
    version,
    long_about = None)]
pub struct Args {
    /// The maximum acceptable marginal false-positive rate.
    #[arg(
        short,
        long = "fp-rate",
        default_value = "1e-6",
        value_name = "FP_RATE"
    )]
    p: f64,
    /// Memory allowance for the Bloom filter, e.g "4GiB". Both binary
    /// (kiB|MiB|GiB) and decimal (kB|MB|GB) formats are understood. Due to
    /// implementation details, a value that is an exact power of 2 (512MiB,
    /// 1GiB, 2GiB etc) gives a modest processing speed advantage (~5%) over
    /// neighbouring values.
    #[arg(short = 'm', long = "mem", default_value = "4GiB")]
    bytes: ByteSize,
    /// Input file. [default: STDIN].
    #[arg(long = "input")]
    input: Option<PathBuf>,
    /// Output file. [default: STDOUT].
    #[arg(long = "output")]
    output: Option<PathBuf>,
    /// Warn instead of error when Bloom filter capacity is exceeded.
    #[arg(long = "allow-overcapacity")]
    allow_overcapacity: bool,
    /// Output metrics file.
    #[arg(long = "metrics", default_value = formatcp!("{PKG_NAME}-metrics.json"))]
    metrics: PathBuf,
    /// Do no work, just print the capacity of the Bloom filter that would be
    /// constructed with the given --fp-rate and --mem values.
    #[arg(long = "show-capacity")]
    show_capacity: bool,
    /// Accept single-ended reads as input [default: paired-end]
    #[arg(long = "single")]
    single: bool,
    /// Unset duplicate flag for any reads that have it set and are no longer
    /// considered duplicate. Only ever required if records have previously
    /// been through a duplicate marking step.
    #[arg(long = "strip-previous")]
    strip_previous: bool,
}

impl Args {
    pub fn reads_per_template(&self) -> u8 {
        match self.single {
            true => 1,
            false => 2,
        }
    }

    #[rustfmt::skip]
    fn show_capacity(&self) {
        println!("specified fp-rate:       {}", self.p);
        println!("specified mem:           {}", self.bytes.to_string_as(true));
        println!("Bloom filter capacity n: {}",
                 BloomFilter::capacity(self.p, self.bytes.as_u64() * 8, DEFAULT_K));
    }
}

pub struct Counts {
    templates: u64,
    templates_unmapped: u64,
    templates_marked_duplicate: u64,
    alignments: u64,
    alignments_marked_duplicate: u64,
}

/// Just logs and exits from unhandled errors
fn main() {
    if let Err(e) = main_() {
        error!("Error: {e}");
        if let Some(source) = e.source() {
            error!("Caused by: {source}");
        }
        std::process::exit(1)
    }
}

/// The main worker
fn main_() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    if args.show_capacity {
        args.show_capacity();
        return Ok(());
    };

    let reader: Reader = match args.input {
        Some(ref path) => Box::new(BufReader::new(
            File::open(path).context(format!("cannot read from {path:?}"))?,
        )),
        None => Box::new(BufReader::new(stdin())),
    };

    let writer: Writer = match args.output {
        Some(ref path) => Box::new(BufWriter::new(
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(path)
                .context(format!("cannot write to {path:?}"))?,
        )),
        None => Box::new(BufWriter::new(stdout())),
    };
    let mut bf = BloomFilter::with_byte_size(args.p, args.bytes.as_u64(), DEFAULT_K);
    info!("BloomFilter capacity: {} items", bf.n());

    let tally = markdups::process_input_stream(reader, writer, &mut bf, &args)?;

    summarize(&args, &bf, &tally)
}

/// Evaluate counts; log summary info and write metrics.
fn summarize(args: &Args, bf: &BloomFilter, tally: &Counts) -> Result<()> {
    let nadded = tally.templates - tally.templates_unmapped - tally.templates_marked_duplicate;
    let fill = 100.0 * nadded as f64 / bf.n() as f64;
    // An implementation of fp on BloomFilter itself using count_estimate for n
    // is very slow for useful values of m.
    let (m, k, n) = (bf.m() as f64, bf.k() as f64, nadded as f64);
    let fp = (1.0 - (-k * n / m).exp()).powf(k);

    let msg_fill = format!("{} items => BloomFilter at {:.2}% capacity", nadded, fill);
    let msg_ok = format!( "Estimated marginal false positive rate: {:.2e}", fp);
    let msg_exc = format!("Target marginal false positive rate {} exceeded", bf.p());

    if fill <= 100.0 {
        info!("{msg_fill}");
        info!("{msg_ok}");
    } else if args.allow_overcapacity {
        warn!("{msg_fill}");
        warn!("{msg_exc}");
    } else {
        error!("{msg_fill}");
        bail!("{msg_exc}");
    }

    let dup_frac = (tally.templates_marked_duplicate as f64)
        / (tally.templates - tally.templates_unmapped) as f64;
    let metrics = json!({
        "ALIGNMENTS": tally.alignments,
        "ALIGNMENTS_MARKED_DUPLICATE": tally.alignments_marked_duplicate,
        "TEMPLATES": tally.templates,
        "TEMPLATES_UNMAPPED": tally.templates_unmapped,
        "TEMPLATES_MARKED_DUPLICATE": tally.templates_marked_duplicate,
        "TEMPLATE_DUPLICATE_FRACTION": dup_frac
    });

    info!("{}", to_string_pretty(&metrics).unwrap());

    OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&args.metrics)
        .context(format!("cannot write to {:?}", args.metrics))?;
    Ok(())
}
