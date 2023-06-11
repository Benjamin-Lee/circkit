use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

pub fn orfs(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Orfs {
            input,
            output,
            min_length,
            start_codons,
            stop_codons,
            include_stop,
            no_wrap,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            // Step 1: Find all stop and start codons by frame
            let start_codons = start_codons.split(',').collect::<Vec<_>>();
            let stop_codons = stop_codons.split(',').collect::<Vec<_>>();

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

                    let (starts, stops) = circkit::orfs::start_stop_codon_indices_by_frame_naive(
                        &std::str::from_utf8(&normalized).unwrap(),
                        &start_codons,
                        &stop_codons,
                    );

                    let mut all_orfs =
                        circkit::orfs::find_orfs_with_indices(normalized.len(), starts, stops);

                    // length filtering
                    all_orfs.retain(|orf| orf.length >= *min_length);

                    // wrap filtering with divisibility by 3
                    if *no_wrap {
                        all_orfs.retain(|orf| {
                            orf.start < orf.stop.unwrap_or(orf.length) && orf.length % 3 == 0
                        });
                    }

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
                            .write_all(
                                orf.seq_with_opts(&record.full_seq(), *include_stop)
                                    .as_bytes(),
                            )
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
