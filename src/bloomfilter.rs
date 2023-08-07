use arrayvec::ArrayVec;
use bit_vec::BitVec;
use bytesize::ByteSize;
use log::{debug, info};
use xxhash_rust::xxh3::xxh3_64_with_seed;

// any two primes "should" do for double hashing
const DEFAULT_SEED1: u64 = 43;
const DEFAULT_SEED2: u64 = 9967;
const U8_RANGE: usize = u8::MAX as usize + 1;

#[derive(Debug)]
pub struct BloomFilter {
    /// maximum acceptable marginal false-positive rate
    p: f64,
    /// storage size in bits
    m: u64,
    /// number of hash functions
    k: u8,
    mpow2: bool,
    mask: u64,
    seed1: u64,
    seed2: u64,
    bits: BitVec,
}

#[allow(dead_code)]
impl BloomFilter {
    /*
     * Associated functions
     */

    /// Return the filter capacity `n` inferred from `p`, `m`, and `k`.
    pub fn capacity(p: f64, m: u64, k: u8) -> u64 {
        let (mf, kf) = (m as f64, k as f64);
        (mf / (-kf / (1.0 - (p.ln() / kf).exp()).ln())).ceil() as u64
    }

    /// Return the memory-optimal bitset size `m` and number of hash functions
    /// `k` for a given number of items to store `n` and target false positive
    /// rate `p`.
    pub fn m_k_min(p: f64, n: u64) -> (u64, u8) {
        let m = ((n as f64 * -p.ln()) / (2.0_f64.ln()).powf(2.0)).ceil() as u64;
        let k = (2.0_f64.ln() * m as f64 / n as f64).ceil() as u8;
        (m, k)
    }

    /// Create a Bloom filter with false-positive rate `p`, bit array size `m`,
    /// and `k` hash functions.
    pub fn new(p: f64, m: u64, k: u8) -> Self {
        let mpow2 = m & (m - 1) == 0;
        debug!(
            "mpow2={0}; {1} will be used for hash addressing",
            mpow2,
            if mpow2 { "bit mask" } else { "modulus" }
        );
        let bf = Self {
            p,
            m,
            k,
            mpow2,
            mask: if mpow2 { m - 1 } else { 0 },
            seed1: DEFAULT_SEED1,
            seed2: DEFAULT_SEED2,
            bits: BitVec::from_elem(m as usize, false),
        };
        info!("BloomFilter initialized with p={0} m={1} k={2}", p, m, k);
        bf
    }

    /// Create a Bloom filter with false-positive rate `p`, bit array size
    /// specified as `nbytes`, and `k` hash functions.
    pub fn with_byte_size(p: f64, nbytes: u64, k: u8) -> Self {
        Self::new(p, 8 * nbytes, k)
    }

    /// Create a Bloom filter with false-positive rate `p` and capacity `n`.
    /// `k` is determined such that bit array `m` size is minimized.
    pub fn with_capacity(p: f64, n: u64) -> Self {
        let (m, k) = Self::m_k_min(p, n);
        Self::new(p, m, k)
    }

    /// Create a Bloom filter with false-positive rate `p`, bit array size
    /// specified as a human-readable byte size string `mem`, and `k` hash
    /// functions.
    pub fn with_mem_spec(p: f64, mem: &str, k: u8) -> Self {
        Self::with_byte_size(p, mem.parse::<ByteSize>().unwrap().as_u64(), k)
    }

    /*
     * public methods
     */

    /// Add the item.
    pub fn add(&mut self, item: &[u8]) {
        // Profiling shows:
        // * using ArrayVec is substantially faster than Vec::with_capacity
        //   because of avoided alloc/free
        // * declaring offsets here and passing it by reference to hash() is
        //   faster than having hash() declare and return it
        let mut offsets = ArrayVec::<_, U8_RANGE>::new();
        self.hash(item, &mut offsets);
        for x in offsets.iter() {
            self.bits.set(*x, true);
        }
    }

    /// Add the item; return false if it was already present otherwise true.
    pub fn add_once(&mut self, item: &[u8]) -> bool {
        let mut offsets = ArrayVec::<_, U8_RANGE>::new();
        self.hash(item, &mut offsets);
        // This is optimized for applications where mostly unique items are
        // added: all() is short-circuiting. For applications where duplicates
        // predominate it may be cheaper to test every bit and set if unset.
        if offsets.iter().all(|x| self.bits.get(*x).unwrap()) {
            false
        } else {
            for x in offsets.iter() {
                self.bits.set(*x, true);
            }
            true
        }
    }

    /// Check if the item is present.
    pub fn check(&self, item: &[u8]) -> bool {
        let mut offsets = ArrayVec::<_, U8_RANGE>::new();
        self.hash(item, &mut offsets);
        offsets.iter().all(|x| self.bits.get(*x).unwrap())
    }

    /// Return estimated number of items stored. The current implementation
    /// is very expensive because it accesses memory `m` times. Probably not
    /// useful for most values of `m`.
    ///
    /// Ref:Swamidass & Baldi (2007) <https://doi.org/10.1021/ci600358f>
    pub fn count_estimate(&self) -> u64 {
        let (m_, k_) = (self.m as f64, self.k as f64);
        let nset = self.bits.iter().filter(|x| *x).count() as f64;
        ((m_ / k_) * -(1.0 - nset / m_).ln()).ceil() as u64
    }

    /// filter capacity `n` inferred from `p`, `m`, and `k`.
    pub fn n(&self) -> u64 {
        Self::capacity(self.p, self.m, self.k)
    }

    /// filter false positive rate `p`
    pub fn p(&self) -> f64 {
        self.p
    }

    /// filter bit array size `m`
    pub fn m(&self) -> u64 {
        self.m
    }

    /// filter hashes `k`
    pub fn k(&self) -> u8 {
        self.k
    }

    /*
     * private methods
     */

    /// Generate hash values and map onto k bitarray offsets.
    ///
    /// `k` linear combinations of just 2 independent hashes ("double hashing")
    /// has the same asymptotic behaviour as `k` independent hashes.
    ///
    /// Ref: Kirsch & Mitzenmacher (2006) <https://doi.org/10.1007/11841036_42>
    fn hash(&self, item: &[u8], offsets: &mut ArrayVec<usize, U8_RANGE>) {
        let mut a = xxh3_64_with_seed(item, self.seed1);
        let mut b = xxh3_64_with_seed(item, self.seed2);
        for i in 0..self.k {
            offsets.push(self.offset(&a));
            a = a.wrapping_add(b);
            b = b.wrapping_add(i as u64);
        }
    }

    /// Map a u64 value to a bitarray offset
    fn offset(&self, h: &u64) -> usize {
        // use bitmask if m is a power of 2, otherwise use modulus
        (if self.mpow2 {
            h & self.mask
        } else {
            h % self.m
        }) as usize
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use pretty_assertions::assert_eq;

    // handy conversion of human-readable byte size to bit size
    fn bit_size(s: &str) -> u64 {
        8 * s.parse::<ByteSize>().unwrap().as_u64()
    }

    /// Capacity is calculated correctly e.g. https://hur.st/bloomfilter
    #[test]
    fn capacity() {
        // p=1e-6, m=8MiB, k=6 => n ~ 1 million items
        assert_eq!(BloomFilter::capacity(1e-6, 67108864, 6), 1_178_438);
    }

    /// `k` and minimal `m` are calculated correctly e.g. https://hur.st/bloomfilter
    #[test]
    fn m_k_min() {
        assert_eq!(BloomFilter::m_k_min(1e-6, 1_000_000), (28_755_176, 20));
        assert_eq!(BloomFilter::m_k_min(1e-7, 10_000_000), (335_477_044, 24));
        assert_eq!(BloomFilter::m_k_min(1e-8, 100_000_000), (3_834_023_351, 27));
        assert_eq!(
            BloomFilter::m_k_min(1e-6, 1_000_000_000),
            (28_755_175_133, 20)
        );
    }

    #[test]
    fn new_mpow2_true() {
        let (p, m, k) = (1e-6, 33_554_432, 10);
        assert!(m == bit_size("4MiB"));
        let bf = BloomFilter::new(p, m, k);
        assert_eq!(bf.p(), p);
        assert_eq!(bf.m(), m);
        assert_eq!(bf.k(), k);
        assert_eq!(bf.mpow2, true);
    }

    #[test]
    fn new_mpow2_false() {
        let (p, m, k) = (1e-6, 32_000_000, 10);
        assert!(m == bit_size("4MB"));
        let bf = BloomFilter::new(p, m, k);
        assert_eq!(bf.p(), p);
        assert_eq!(bf.m(), m);
        assert_eq!(bf.k(), k);
        assert_eq!(bf.mpow2, false);
    }

    #[test]
    fn with_byte_size() {
        let (p, m, k) = (1e-6, 33_554_432, 10);
        assert_eq!(BloomFilter::with_byte_size(p, m / 8, k).m(), m);
    }

    #[test]
    fn with_capacity() {
        let (p, n) = (1e-6, 1_000_000);
        let bf = BloomFilter::with_capacity(p, n);
        assert_eq!(bf.m(), 28_755_176);
        assert_eq!(bf.k(), 20);
    }

    #[test]
    fn with_mem_spec() {
        let (p, m, k) = (1e-6, 33_554_432, 10);
        assert_eq!(BloomFilter::with_mem_spec(p, "4MiB", k).m(), m);
    }

    /// add_once() returns true when the item was absent.
    #[test]
    fn add_once_missing() {
        let mut bf = BloomFilter::new(0.001, 1000, 2);
        let key = "something".as_bytes();
        assert_eq!(bf.add_once(key), true);
    }

    /// add_once() returns false when the item was present.
    #[test]
    fn add_once_existing() {
        let mut bf = BloomFilter::new(0.001, 1000, 2);
        let key = "something".as_bytes();
        bf.add(key);
        assert_eq!(bf.add_once(key), false);
    }

    /// check() returns false when the item is absent.
    #[test]
    fn check_missing() {
        let bf = BloomFilter::new(0.001, 1000, 2);
        let key = "something".as_bytes();
        assert_eq!(bf.check(key), false);
    }

    /// check() returns true when the item is present.
    #[test]
    fn check_existing() {
        let mut bf = BloomFilter::new(0.001, 1000, 2);
        let key = "something".as_bytes();
        bf.add(key);
        assert_eq!(bf.check(key), true);
    }

    /// Estimate for number of stored items is within 0.001 of true value.
    #[test]
    fn count_estimate() {
        let mut bf = BloomFilter::with_mem_spec(1e-6, "8MiB", 6);
        let nitems = 1_000_000;
        for i in 0..nitems {
            bf.add(&i.to_string().as_bytes());
        }
        assert!((1.0 - (bf.count_estimate() as f64 / nitems as f64)).abs() < 0.001);
    }

    /// False negative rate should be exactly zero.
    #[test]
    fn fnr_eq_0() {
        let mut bf = BloomFilter::with_mem_spec(1e-6, "8MiB", 6);
        let mut not_present = 0;
        let n = 1_000_000;
        for i in 0..n {
            let value = i.to_string();
            bf.add(value.as_bytes());
            if !bf.check(value.as_bytes()) {
                not_present += 1
            }
        }
        assert_eq!(not_present, 0);
    }

    /// At capacity, the marginal false-positive rate is less than 1.1*p.
    #[test]
    fn fpr_bound() {
        for exp in -7..=-1 {
            let p = 10_f64.powi(exp);
            let mut bf = BloomFilter::with_mem_spec(p, "8MiB", 6);
            let n = bf.n() as i64;
            let mut fp = 0;
            for i in 0..n {
                // true positives are 0 and the +ve integers as native-endian bytes
                bf.add(&i.to_ne_bytes());
            }
            for i in 1..=n {
                // true negatives are the -ve integers as native-endian bytes
                if bf.check(&(-i).to_ne_bytes()) {
                    fp += 1;
                }
            }
            let fpr = fp as f64 / n as f64;
            // allow 10% wiggle room: p is not a strict bound
            assert!(fpr / (p) <= 1.1, "p: {p}, n: {n}, fpr: {fpr}");
        }
    }
}
