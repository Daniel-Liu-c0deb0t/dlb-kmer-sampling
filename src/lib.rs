use rustc_hash::{FxHashMap, FxHashSet};

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

const MIN_S: usize = 5;

impl Seeder {
    pub fn new_optimal_closed_syncmer(k: usize, w: usize) -> Self {
        assert!(k < 64);
        assert!(k % 2 != 0);

        let mut t = 0;

        while k >= MIN_S + (t + 1) * w {
            t += 1;
        }

        assert!(t > 0);

        let s = k - t * w;
        let mut mask = 0u64;

        for i in 0..=t {
            mask |= 1 << (i * w);
        }

        Self { k, s, t, mask }
    }

    pub fn new_closed_syncmer(k: usize, w: usize) -> Self {
        assert!(k < 64);
        assert!(k % 2 != 0);

        let s = k - w;
        let mut mask = 0u64;

        for i in 0..=1 {
            mask |= 1 << (i * w);
        }

        Self { k, s, t: 1, mask }
    }

    pub fn new_distance_lower_bound(k: usize, reciprocal_density: usize) -> Self {
        assert!(k < 64);
        assert!(k % 2 != 0);

        let mut t = 0;
        let mut best_t = 0;
        let mut best_stride = 0;

        while k + 1 >= MIN_S + (t + 1) * reciprocal_density {
            t += 1;

            let curr_stride = (reciprocal_density - 1) * t / (t + 1);

            if t == 1 || curr_stride > best_stride {
                best_t = t;
                best_stride = curr_stride;
            }
        }

        let t = best_t;

        assert!(t > 0);

        let s = k - t * reciprocal_density + 1;
        let mut mask = 0u64;
        let stride = (k - s + 1 - t) / (t + 1);

        for i in 0..t {
            mask |= 1 << ((i + 1) * stride + i);
        }

        Self { k, s, t, mask }
    }

    pub fn distance_lower_bound(&self) -> usize {
        (self.k - self.s + 1 - self.t) / (self.t + 1) + 1
    }

    pub fn new_start(k: usize, reciprocal_density: usize) -> Self {
        Self {
            k,
            s: k - reciprocal_density + 1,
            t: 0,
            mask: 0b1,
        }
    }

    pub fn new_mid(k: usize, reciprocal_density: usize) -> Self {
        Self {
            k,
            s: k - reciprocal_density + 1,
            t: 1,
            mask: 1 << ((reciprocal_density - 1) / 2),
        }
    }

    pub fn get_seeds(&self, seq: &[u8], res: &mut Vec<(u128, usize)>) {
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

pub struct ModMinimizer {
    pub k: usize,
    pub w: usize,
    pub s: usize,
    pub l: usize,
}

const MIN_S_MOD: usize = 4;

impl ModMinimizer {
    pub fn new_original(k: usize, w: usize) -> Self {
        let l = k + w - 1;
        let s = MIN_S_MOD + (k - MIN_S_MOD) % w;

        Self { k, w, s, l }
    }

    pub fn new_dw(k: usize, w: usize) -> Self {
        let l = k + w - 1;

        let density = |k, w| {
            let t = (k + w - MIN_S_MOD) / w;
            ((t + 1) as f64) / ((t * w + 1) as f64)
        };

        let mut w2 = w;
        let mut best_w = 0;
        let mut best_density = std::f64::MAX;

        while w2 > 0 && w2 >= w / 2 {
            let curr_density = density(k, w2);

            if curr_density < best_density {
                best_density = curr_density;
                best_w = w2;
            }

            w2 -= 1;
        }

        let best_k = l + 1 - best_w;
        let s = MIN_S_MOD + (best_k - MIN_S_MOD) % best_w;

        Self { k, w: best_w, s, l }
    }

    pub fn get_seeds(&self, seq: &[u8], res: &mut Vec<(Vec<u8>, usize)>) {
        if seq.len() < self.l {
            return;
        }

        for (i, win) in seq.windows(self.l).enumerate() {
            let mut s_min = std::u64::MAX;
            let mut s_idx = 0;

            for j in 0..=(self.l - self.s) {
                let s_mer = &win[j..j + self.s];
                let s_hash = hash_slice(0, s_mer);

                if s_hash <= s_min {
                    s_min = s_hash;
                    s_idx = j;
                }
            }

            let k_idx = s_idx % self.w;
            let kmer = win[k_idx..k_idx + self.k].to_vec();

            if let Some(last) = res.last() {
                let (last_kmer, last_idx) = last;

                if &kmer != last_kmer || i + k_idx != *last_idx {
                    res.push((kmer, i + k_idx));
                }
            } else {
                res.push((kmer, i + k_idx));
            }
        }
    }
}

pub struct SubsampledMinimizer {
    pub l: usize,
    pub k: usize,
    pub w: usize,
    pub s: usize,
}

const MIN_S_SUB: usize = 3;

impl SubsampledMinimizer {
    pub fn new(k: usize, w: usize) -> Self {
        let l = k + w - 1;
        let s = (k + 1).saturating_sub(w).max(MIN_S_SUB);
        Self { l, k, w, s }
    }

    pub fn get_seeds(&self, seq: &[u8], res: &mut Vec<(Vec<u8>, usize)>) {
        if seq.len() < self.l {
            return;
        }

        for (i, win) in seq.windows(self.l).enumerate() {
            let mut k_min = (10, std::u64::MAX);
            let mut k_idx = 0;

            for j in 0..self.w {
                let k_mer = &win[j..j + self.k];
                let mut s_min = std::u64::MAX;
                let mut s_idx = 0;

                for k in 0..=(self.k - self.s) {
                    let s_mer = &k_mer[k..k + self.s];
                    let s_hash = hash_slice(0, s_mer);

                    if s_hash <= s_min {
                        s_min = s_hash;
                        s_idx = k;
                    }
                }

                let mut k_hash = k_min;

                if s_idx == (self.k - self.s + 1) / 2 {
                    k_hash = (0, hash_slice(0, k_mer));
                }

                if s_idx == 0 || s_idx == self.k - self.s {
                    k_hash = (1, hash_slice(0, k_mer));
                }

                if k_hash < k_min {
                    k_min = k_hash;
                    k_idx = j;
                }
            }

            let kmer = win[k_idx..k_idx + self.k].to_vec();

            if let Some(last) = res.last() {
                let (last_kmer, last_idx) = last;

                if &kmer != last_kmer || i + k_idx != *last_idx {
                    res.push((kmer, i + k_idx));
                }
            } else {
                res.push((kmer, i + k_idx));
            }
        }
    }
}

fn hash(seed: u64, kmer: u128) -> u64 {
    let mut h = WyHash::with_seed(seed);
    kmer.hash(&mut h);
    h.finish()
}

fn hash_slice(seed: u64, kmer: &[u8]) -> u64 {
    let mut h = WyHash::with_seed(seed);
    kmer.hash(&mut h);
    h.finish()
}

pub fn density<T>(seeds: &[(T, usize)], seq_len: usize, k: usize) -> f64 {
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

pub fn conservation(
    seeds: &[(u128, usize)],
    seeds_mut: &[(u128, usize)],
    seq_len: usize,
    k: usize,
) -> f64 {
    let mut res = FxHashSet::default();
    let a = seeds.iter().clone().collect::<FxHashSet<_>>();
    let b = seeds_mut.iter().clone().collect::<FxHashSet<_>>();

    for &&(_, i) in a.intersection(&b) {
        for j in i..i + k {
            res.insert(j);
        }
    }

    (res.len() as f64) / (seq_len as f64)
}

pub fn l2(seeds: &[(u128, usize)], seeds_mut: &[(u128, usize)], seq_len: usize, k: usize) -> f64 {
    let a = seeds.iter().clone().collect::<FxHashSet<_>>();
    let b = seeds_mut.iter().clone().collect::<FxHashSet<_>>();
    let mut c = a.intersection(&b).map(|&&s| s).collect::<Vec<_>>();
    c.sort_unstable_by_key(|(_, i)| *i);
    let mut res = 0.0f64;

    for w in c.windows(2) {
        let &[(_, i), (_, j)] = w else { unreachable!() };
        res += (j - i).saturating_sub(k).pow(2) as f64;
    }

    (res / (seq_len as f64)).sqrt()
}
