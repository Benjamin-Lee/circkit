use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use aho_corasick::AhoCorasick;
use seq_io::{
    fasta::{Record, RefRecord},
    parallel::parallel_fasta,
};

pub fn orfs(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Orfs {
            input,
            output,
            min_length,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            // Step 1: Find all stop and start codons by frame
            let start_codons = ["ATG"];
            let stop_codons = ["TAA", "TAG", "TGA"];
            let patterns = start_codons
                .iter()
                .chain(stop_codons.iter())
                .map(|s| s.as_bytes())
                .collect::<Vec<_>>();
            let ac = AhoCorasick::new(patterns).unwrap();

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, orfs| {
                    // runs in worker
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    // get start/stop codon indexing method from INDEX_METHOD env var
                    let index_method = std::env::var("INDEX_METHOD").unwrap_or("aho".to_string());

                    let (starts, stops) = match index_method.as_str() {
                        "naive" => circkit::orfs::start_stop_codon_indices_by_frame_naive(
                            &std::str::from_utf8(&normalized).unwrap(),
                            &start_codons,
                            &stop_codons,
                        ),
                        "iter" => circkit::orfs::start_stop_codon_indices_by_frame_iter(
                            &std::str::from_utf8(&normalized).unwrap(),
                            &start_codons,
                            &stop_codons,
                        ),
                        "aho" => circkit::orfs::start_stop_codon_indices_by_frame_aho_corasick(
                            &std::str::from_utf8(&normalized).unwrap(),
                            &start_codons,
                            &stop_codons,
                            &ac,
                        ),
                        _ => panic!("INDEX_METHOD env var must be either 'iter' or 'aho'"),
                    };

                    let mut all_orfs =
                        circkit::orfs::find_orfs_with_indices(normalized.len(), starts, stops);

                    // length filtering
                    all_orfs.retain(|orf| orf.length >= *min_length);

                    *orfs = circkit::orfs::longest_orfs(&mut all_orfs);
                },
                |record, orfs| {
                    for orf in orfs {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.head()).unwrap();
                        writer.write_all(b" ORF").unwrap();
                        writer.write_all(orf.start.to_string().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer
                            .write_all(orf.seq(&record.full_seq()).as_bytes())
                            .unwrap();
                        writer.write_all(b"\n").unwrap();
                    }

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
        }
        _ => panic!("input command is not for orfs"),
    }
    Ok(())
}
