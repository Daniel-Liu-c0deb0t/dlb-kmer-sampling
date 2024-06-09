use simulate_seqs::*;

use dlb_kmer_sampling::*;

use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let reciprocal_density = 11;
    let w = 11;
    let ks = [21, 31, 41, 51, 61];
    let mut_rates = [0.01, 0.05, 0.1, 0.15];

    let mut distances_file = BufWriter::new(File::create("distances.csv").unwrap());
    let mut distances_closed_file = BufWriter::new(File::create("distances_closed.csv").unwrap());
    let mut metrics_file = BufWriter::new(File::create("metrics.csv").unwrap());
    let mut metrics_closed_file = BufWriter::new(File::create("metrics_closed.csv").unwrap());

    let mut rng = StdRng::seed_from_u64(0);
    let alpha = [b'A', b'C', b'G', b'T'];
    let seq = rand_str(1000000, &alpha, &mut rng);
    let seqs_mut = mut_rates
        .iter()
        .map(|&r| rand_mutate_subs(&seq, ((seq.len() as f64) * r) as usize, &alpha, &mut rng))
        .collect::<Vec<_>>();

    writeln!(&mut distances_file, "algorithm,k,distance,count").unwrap();
    writeln!(
        &mut distances_closed_file,
        "algorithm,k,density,distance,count"
    )
    .unwrap();
    writeln!(&mut metrics_file, "algorithm,k,mut_rate,metric,value").unwrap();
    writeln!(
        &mut metrics_closed_file,
        "algorithm,k,mut_rate,metric,value"
    )
    .unwrap();

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
            writeln!(
                &mut distances_file,
                "Open syncmer (Edgar 2021),{k},{dist},{count}"
            )
            .unwrap();
        }

        for (&r, seq_mut) in mut_rates.iter().zip(&seqs_mut) {
            let mut seeds_mut = Vec::new();
            seeder.get_seeds(&seq_mut, &mut seeds_mut);
            let conservation = conservation(&seeds, &seeds_mut, seq.len(), k);
            let l2 = l2(&seeds, &seeds_mut, seq.len(), k);
            writeln!(
                &mut metrics_file,
                "Open syncmer (Edgar 2021),{k},{r},Conservation,{conservation}"
            )
            .unwrap();
            writeln!(
                &mut metrics_file,
                "Open syncmer (Edgar 2021),{k},{r},L2,{l2}"
            )
            .unwrap();
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
            writeln!(
                &mut distances_file,
                "Open syncmer (Shaw & Yu 2022),{k},{dist},{count}"
            )
            .unwrap();
        }

        for (&r, seq_mut) in mut_rates.iter().zip(&seqs_mut) {
            let mut seeds_mut = Vec::new();
            seeder.get_seeds(&seq_mut, &mut seeds_mut);
            let conservation = conservation(&seeds, &seeds_mut, seq.len(), k);
            let l2 = l2(&seeds, &seeds_mut, seq.len(), k);
            writeln!(
                &mut metrics_file,
                "Open syncmer (Shaw & Yu 2022),{k},{r},Conservation,{conservation}"
            )
            .unwrap();
            writeln!(
                &mut metrics_file,
                "Open syncmer (Shaw & Yu 2022),{k},{r},L2,{l2}"
            )
            .unwrap();
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
            writeln!(&mut distances_file, "DLB (ours),{k},{dist},{count}").unwrap();
        }

        for (&r, seq_mut) in mut_rates.iter().zip(&seqs_mut) {
            let mut seeds_mut = Vec::new();
            seeder.get_seeds(&seq_mut, &mut seeds_mut);
            let conservation = conservation(&seeds, &seeds_mut, seq.len(), k);
            let l2 = l2(&seeds, &seeds_mut, seq.len(), k);
            writeln!(
                &mut metrics_file,
                "DLB (ours),{k},{r},Conservation,{conservation}"
            )
            .unwrap();
            writeln!(&mut metrics_file, "DLB (ours),{k},{r},L2,{l2}").unwrap();
        }
    }

    for &k in &ks {
        let seeder = Seeder::new_closed_syncmer(k, w);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        eprintln!("Closed syncmer (Edgar 2021)");
        eprintln!("k: {k}");
        eprintln!("s: {}", seeder.s);
        eprintln!("density: {density}");
        eprintln!("distances: {distances:?}");
        eprintln!();

        for (dist, count) in distances {
            writeln!(
                &mut distances_closed_file,
                "Closed syncmer (Edgar 2021),{k},{density},{dist},{count}"
            )
            .unwrap();
        }

        for (&r, seq_mut) in mut_rates.iter().zip(&seqs_mut) {
            let mut seeds_mut = Vec::new();
            seeder.get_seeds(&seq_mut, &mut seeds_mut);
            let conservation = conservation(&seeds, &seeds_mut, seq.len(), k);
            let l2 = l2(&seeds, &seeds_mut, seq.len(), k);
            writeln!(
                &mut metrics_closed_file,
                "Closed syncmer (Edgar 2021),{k},{r},Conservation,{conservation}"
            )
            .unwrap();
            writeln!(
                &mut metrics_closed_file,
                "Closed syncmer (Edgar 2021),{k},{r},L2,{l2}"
            )
            .unwrap();
        }
    }

    for &k in &ks {
        let seeder = Seeder::new_optimal_closed_syncmer(k, w);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);
        let distances = distances(&seeds);

        eprintln!("Optimal closed syncmer (ours)");
        eprintln!("k: {k}");
        eprintln!("s: {}", seeder.s);
        eprintln!("t: {}", seeder.t);
        eprintln!("density: {density}");
        eprintln!("distances: {distances:?}");
        eprintln!();

        for (dist, count) in distances {
            writeln!(
                &mut distances_closed_file,
                "Optimal closed syncmer (ours),{k},{density},{dist},{count}"
            )
            .unwrap();
        }

        for (&r, seq_mut) in mut_rates.iter().zip(&seqs_mut) {
            let mut seeds_mut = Vec::new();
            seeder.get_seeds(&seq_mut, &mut seeds_mut);
            let conservation = conservation(&seeds, &seeds_mut, seq.len(), k);
            let l2 = l2(&seeds, &seeds_mut, seq.len(), k);
            writeln!(
                &mut metrics_closed_file,
                "Optimal closed syncmer (ours),{k},{r},Conservation,{conservation}"
            )
            .unwrap();
            writeln!(
                &mut metrics_closed_file,
                "Optimal closed syncmer (ours),{k},{r},L2,{l2}"
            )
            .unwrap();
        }
    }
}
