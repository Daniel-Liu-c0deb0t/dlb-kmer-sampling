use rustc_hash::FxHashMap;

use wyhash::WyHash;

use std::hash::{Hash, Hasher};

static BYTE2BITS: [u8; 256] = {
    let mut l = [0u8; 256];
    l[b'A' as usize] = 0b00;
    l[b'a' as usize] = 0b00;
    l[b'C' as usize] = 0b01;
    l[b'c' as usize] = 0b01;
    l[b'G' as usize] = 0b10;
    l[b'g' as usize] = 0b10;
    l[b'T' as usize] = 0b11;
    l[b't' as usize] = 0b11;
    l
};

pub struct Seeder {
    pub k: usize,
    pub s: usize,
    pub t: usize,
    pub mask: u64,
}

const MIN_S: usize = 7;

impl Seeder {
    pub fn new(k: usize, reciprocal_density: usize) -> Self {
        assert!(k < 64);
        assert!(k % 2 != 0);

        let mut t = 0;

        while k + 1 >= MIN_S + (t + 1) * reciprocal_density
            && (t + 1) * (reciprocal_density - 1) >= t + 2 {
            t += 1;
        }

        assert!(t > 0);

        let s = k - t * reciprocal_density + 1;
        let mut mask = 0u64;
        let stride = (k - s + 1 - t) / (t + 1);

        for i in 0..t {
            mask |= 1 << ((i + 1) * stride + i);
        }

        Self { k, s, t, mask }
    }

    pub fn get_seeds(
        &self,
        seq: &[u8],
        res: &mut Vec<(u128, usize)>,
    ) {
        if seq.len() < self.k {
            return;
        }

        let s_mask = (1u128 << (self.s * 2)) - 1;
        let k_mask = (1u128 << (self.k * 2)) - 1;
        let mut kmer = 0u128;

        for &next in &seq[..self.k - 1] {
            unsafe {
                kmer = (kmer << 2) | (*BYTE2BITS.as_ptr().add(next as usize) as u128);
            }
        }

        for (i, win) in seq.windows(self.k).enumerate() {
            let next = *win.last().unwrap();

            unsafe {
                kmer = ((kmer << 2) & k_mask) | (*BYTE2BITS.as_ptr().add(next as usize) as u128);
            }

            let mut s_min = std::u64::MAX;
            let mut s_idx = 0;

            for j in 0..=(self.k - self.s) {
                let s_mer = (kmer >> ((self.k - self.s - j) * 2)) & s_mask;
                let s_hash = hash(0, s_mer);

                if s_hash <= s_min {
                    s_min = s_hash;
                    s_idx = j;
                }
            }

            if (self.mask >> s_idx) & 0b1 > 0 {
                res.push((kmer, i));
            }
        }
    }
}

fn hash(seed: u64, kmer: u128) -> u64 {
    let mut h = WyHash::with_seed(seed);
    kmer.hash(&mut h);
    h.finish()
}

pub fn density(seeds: &[(u128, usize)], seq_len: usize, k: usize) -> f64 {
    if seq_len < k {
        return 0.0;
    }

    (seeds.len() as f64) / ((seq_len - k + 1) as f64)
}

pub fn distances(seeds: &[(u128, usize)]) -> Vec<(usize, usize)> {
    let mut seeds = seeds.to_owned();
    seeds.sort_unstable_by_key(|(_, i)| *i);
    let mut hist = FxHashMap::default();
    let mut prev = seeds[0].1;

    for &(_, i) in &seeds[1..] {
        *hist.entry(i - prev).or_insert(0usize) += 1;
        prev = i;
    }

    let mut res = hist.into_iter().collect::<Vec<_>>();
    res.sort_unstable();
    res
}
