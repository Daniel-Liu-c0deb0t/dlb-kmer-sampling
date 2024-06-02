use simulate_seqs::*;

use conserved_kmer_sampling::*;

fn main() {
    let reciprocal_density = 11;

    let mut rng = StdRng::seed_from_u64(0);
    let alpha = [b'A', b'C', b'G', b'T'];
    let seq = rand_str(100000, &alpha, &mut rng);

    for k in [21, 31, 41, 51, 61] {
        let seeder = Seeder::new(k, reciprocal_density);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        println!("k: {k}");
        println!("t: {}", seeder.t);
        println!("Density: {density}");
        println!("Distances: {distances:?}");
        println!();
    }
}
