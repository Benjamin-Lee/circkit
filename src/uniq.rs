use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use rustc_hash::{FxHashSet, FxHasher};
use seq_io::{fasta::Record, parallel::parallel_fasta};
use std::{
    hash::{Hash, Hasher},
    sync::Mutex,
};

pub fn uniq(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Uniq {
            input,
            output,
            normalize,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            let seen = Mutex::new(FxHashSet::default());

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, seq| {
                    // runs in worker
                    let normalized =
                        circkit::normalize(&record.full_seq(), circkit::normalize::Alphabet::Dna);

                    // hash the normalized sequence
                    let mut hasher = FxHasher::default();
                    normalized.hash(&mut hasher);
                    let hash = hasher.finish();

                    let mut x = seen.lock().unwrap();
                    if !x.insert(hash) {
                        *seq = vec![];
                    } else {
                        // if the user requested, write the normalized sequence
                        if *normalize {
                            *seq = normalized
                        } else {
                            *seq = record.seq().to_vec()
                        }
                    }
                    drop(x);
                },
                |record, seq| {
                    // runs in main thread
                    if !seq.is_empty() {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.id().unwrap().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer.write_all(seq).unwrap();
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
        _ => panic!("input command is not for uniq"),
    }
    Ok(())
}
