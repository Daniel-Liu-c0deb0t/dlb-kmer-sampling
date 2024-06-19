use simulate_seqs::*;

use dlb_kmer_sampling::*;

use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let w = 24;
    let ks = (5..200).collect::<Vec<_>>();

    let mut density_file = BufWriter::new(File::create("density_mod.csv").unwrap());

    let mut rng = StdRng::seed_from_u64(0);
    let alpha = [b'A', b'C', b'G', b'T'];
    let seq = rand_str(1000000, &alpha, &mut rng);

    writeln!(&mut density_file, "algorithm,k,w,density").unwrap();

    for &k in &ks {
        let seeder = ModMinimizer::new_original(k, w);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);

        writeln!(
            &mut density_file,
            "Mod minimizer,{k},{},{density}",
            seeder.w,
        )
        .unwrap();
    }

    for &k in &ks {
        let seeder = ModMinimizer::new_dw(k, w);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);

        writeln!(
            &mut density_file,
            "Mod minimizer DW,{k},{},{density}",
            seeder.w,
        )
        .unwrap();
    }

    for &k in &ks {
        let seeder = SubsampledMinimizer::new(k, w);
        let mut seeds = Vec::new();
        seeder.get_seeds(&seq, &mut seeds);

        let density = density(&seeds, seq.len(), k);

        writeln!(
            &mut density_file,
            "OpenClosed,{k},{},{density}",
            seeder.w,
        )
        .unwrap();
    }
}
