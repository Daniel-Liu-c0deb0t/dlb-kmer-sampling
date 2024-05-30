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
    k: usize,
    m1: usize,
    m2: usize,
}

impl Seeder {
    pub fn new(k: usize, m1: usize, m2: usize) -> Self {
        assert!(m2 <= m1 && m1 <= k && k < 32);
        assert!(k % 2 != 0 && m1 % 2 != 0 && m2 % 2 != 0);
        assert!(m1 <= ((k - m2 + 1) / 2 + m2));
        Self { k, m1, m2 }
    }

    pub fn get_seeds(
        &self,
        seq: &[u8],
        res: &mut Vec<(u64, usize)>,
    ) {
        if seq.len() < self.k {
            return;
        }

        let m1_mask = (1u64 << (self.m1 * 2)) - 1;
        let m2_mask = (1u64 << (self.m2 * 2)) - 1;
        let k_mask = (1u64 << (self.k * 2)) - 1;
        let mid = (self.k - self.m2 + 1) / 2;
        let mut kmer = 0u64;

        for &next in &seq[..self.k - 1] {
            unsafe {
                kmer = (kmer << 2) | (*BYTE2BITS.as_ptr().add(next as usize) as u64);
            }
        }

        for (i, win) in seq.windows(self.k).enumerate() {
            let next = *win.last().unwrap();

            unsafe {
                kmer = ((kmer << 2) & k_mask) | (*BYTE2BITS.as_ptr().add(next as usize) as u64);
            }

            let mut m1_min = (std::u64::MAX, std::u64::MAX);
            let mut m1_idx = 0;

            for j in 0..=(self.k - self.m1) {
                let m1_mer = (kmer >> ((self.k - self.m1 - j) * 2)) & m1_mask;
                let mut m2_mer_hash = std::u64::MAX;

                for l in 0..=(self.m1 - self.m2) {
                    let m2_mer = (m1_mer >> ((self.m1 - self.m2 - l) * 2)) & m2_mask;
                    m2_mer_hash = m2_mer_hash.min(hash(0, m2_mer));
                }

                let m1_hash = (m2_mer_hash, hash(1, m1_mer));

                if m1_hash <= m1_min {
                    m1_min = m1_hash;
                    m1_idx = j;
                }
            }

            if m1_idx == mid {
                res.push((kmer, i));
            }
        }
    }
}

fn hash(seed: u64, kmer: u64) -> u64 {
    let mut h = WyHash::with_seed(seed);
    kmer.hash(&mut h);
    h.finish()
}

pub fn density(seeds: &[(u64, usize)], seq_len: usize, k: usize) -> f64 {
    if seq_len < k {
        return 0.0;
    }

    (seeds.len() as f64) / ((seq_len - k + 1) as f64)
}

pub fn distances(seeds: &[(u64, usize)]) -> Vec<(usize, usize)> {
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
