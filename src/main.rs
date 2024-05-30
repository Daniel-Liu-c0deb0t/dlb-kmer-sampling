use simulate_seqs::*;

use conserved_kmer_sampling::*;

fn main() {
    let k = 31;
    let m1 = 21;
    let m2 = 11;

    let mut rng = StdRng::seed_from_u64(0);
    let alpha = [b'A', b'C', b'G', b'T'];
    let seq = rand_str(10000, &alpha, &mut rng);

    let seeder = Seeder::new(k, m1, m2);
    let mut seeds = Vec::new();
    seeder.get_seeds(&seq, &mut seeds);

    let density = density(&seeds, seq.len(), k);
    let distances = distances(&seeds);

    println!("Density: {density}");
    println!("Distances: {distances:?}");
}
