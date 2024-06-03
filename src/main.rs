use simulate_seqs::*;

use dlb_kmer_sampling::*;

fn main() {
    let reciprocal_density = 11;
    let ks = [21, 31, 41, 51, 61];

    let mut rng = StdRng::seed_from_u64(0);
    let alpha = [b'A', b'C', b'G', b'T'];
    let seq = rand_str(100000, &alpha, &mut rng);

    println!("algorithm,k,distance,count");

    for &k in &ks {
        let seeder = Seeder::new_start(k, reciprocal_density);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        eprintln!("Open syncmer (Edgar 2021)");
        eprintln!("k: {k}");
        eprintln!("s: {}", seeder.s);
        eprintln!("density: {density}");
        eprintln!("distances: {distances:?}");
        eprintln!();

        for (dist, count) in distances {
            println!("Open syncmer (Edgar 2021),{k},{dist},{count}");
        }
    }

    for &k in &ks {
        let seeder = Seeder::new_mid(k, reciprocal_density);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        eprintln!("Open syncmer (Shaw & Yu 2022)");
        eprintln!("k: {k}");
        eprintln!("s: {}", seeder.s);
        eprintln!("density: {density}");
        eprintln!("distances: {distances:?}");
        eprintln!();

        for (dist, count) in distances {
            println!("Open syncmer (Shaw & Yu 2022),{k},{dist},{count}");
        }
    }

    for &k in &ks {
        let seeder = Seeder::new_distance_lower_bound(k, reciprocal_density);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        eprintln!("DLB (ours)");
        eprintln!("k: {k}");
        eprintln!("s: {}", seeder.s);
        eprintln!("t: {}", seeder.t);
        eprintln!("distance lower bound: {}", seeder.distance_lower_bound());
        eprintln!("density: {density}");
        eprintln!("distances: {distances:?}");
        eprintln!();

        for (dist, count) in distances {
            println!("DLB (ours),{k},{dist},{count}");
        }
    }
}
