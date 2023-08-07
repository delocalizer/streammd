use anyhow::Result;
use once_cell::sync::Lazy;
use regex::Regex;

use super::*;

const SAM_DELIMITER: char = '\t';

/// A highly specialized representation of a SAM record, for the purposes of
/// duplicate marking — probably not generally useful for other things...
struct SamRecord {
    buffer: String,
    qname0: usize,
    qname1: usize,
    flag0: usize,
    flag1: usize,
    flag: u16,
    rname0: usize,
    rname1: usize,
    pos: i32,
    cigar0: usize,
    cigar1: usize,
    pg0: usize,
    pg1: usize,
}

impl SamRecord {
    /// Assuming valid boundaries `start` and `stop` for the current field in a
    /// SAM record represented by `buffer`, update them to the next field. The
    /// next field must exist; a `panic` ensues if no delimiters remain. When
    /// `start = stop = 0` the updated boundaries are those of the first field.
    fn advance(buffer: &String, start: &mut usize, stop: &mut usize) {
        if *stop != 0 {
            *start = *stop + 1;
        }
        // &str.find(pattern) benches faster than regex when pattern is char
        *stop = *start
            + &buffer[*start..]
                .find(SAM_DELIMITER)
                .expect("find(SAM_DELIMITER) is None");
    }

    fn from_string(buffer: String) -> Result<Self> {
        static RE_PGTAG: Lazy<Regex> = Lazy::new(|| Regex::new("\tPG:Z:").unwrap());

        // we need field positions to do efficient updates to the buffer so a
        // simple split on the delimiter won't do the job. str.match_indices
        // works but is slower than an explicit advance we do here
        let mut start = 0;
        let mut stop = 0;
        let end = buffer.len();

        // QNAME
        Self::advance(&buffer, &mut start, &mut stop);
        let (qname0, qname1) = (start, stop);

        // FLAG
        Self::advance(&buffer, &mut start, &mut stop);
        let (flag0, flag1) = (start, stop);
        let flag = buffer[start..stop]
            .parse::<u16>()
            .context(format!("Bad FLAG in SAM record: {buffer}"))?;

        // RNAME
        Self::advance(&buffer, &mut start, &mut stop);
        let (rname0, rname1) = (start, stop);

        // POS
        Self::advance(&buffer, &mut start, &mut stop);
        let pos = buffer[start..stop]
            .parse::<i32>()
            .context(format!("Bad POS in SAM record: '{buffer}'"))?;

        // CIGAR
        Self::advance(&buffer, &mut start, &mut stop);
        Self::advance(&buffer, &mut start, &mut stop);
        let (cigar0, cigar1) = (start, stop);

        // PG
        let (pg0, pg1) = match RE_PGTAG
            .find(&buffer[stop + 1..])
            .map(|mtch| stop + 1 + mtch.start())
        {
            Some(pgidx) => (
                pgidx,
                buffer[pgidx + 1..]
                    .find(SAM_DELIMITER)
                    .map(|idx| pgidx + 1 + idx)
                    .unwrap_or(end),
            ),
            _ => (end, end),
        };

        Ok(Self {
            buffer,
            qname0,
            qname1,
            flag0,
            flag1,
            flag,
            rname0,
            rname1,
            pos,
            cigar0,
            cigar1,
            pg0,
            pg1,
        })
    }

    fn qname(&self) -> &str {
        &self.buffer[self.qname0..self.qname1]
    }

    fn rname(&self) -> &str {
        &self.buffer[self.rname0..self.rname1]
    }

    fn pg(&self) -> &str {
        &self.buffer[self.pg0..self.pg1]
    }
}

/// Process the input stream. Header lines are written directly to the output
/// stream; reads are dispatched in qname groups for further processing.
pub fn process_input_stream(
    reader: Reader,
    mut writer: Writer,
    bf: &mut BloomFilter,
    args: &Args,
) -> Result<Counts> {
    let mut qname_prev = String::with_capacity(128);
    for line in reader.lines() {
        let record = SamRecord::from_string(line?)?;
        bf.add_once(record.qname().as_bytes());
//        writer.write(record.pg().as_bytes())?;
//        qname_prev = record.qname().to_string();
//        writer.write(b"\n")?;
    }
    writer.flush()?;
    Ok(Counts {
        templates: 1_000_000_000,
        templates_unmapped: 0,
        templates_marked_duplicate: 33,
        alignments: 0,
        alignments_marked_duplicate: 0,
    })
}
