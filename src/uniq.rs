use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer, table_path_to_writer},
};
use nohash_hasher::BuildNoHashHasher;
use seq_io::{fasta::Record, parallel::parallel_fasta};
use std::collections::HashMap;

#[derive(serde::Serialize)]
struct Row<'a> {
    id: &'a str,
    duplicate_id: &'a str,
}

pub fn uniq(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Uniq {
            input,
            output,
            canonicalize,
            table,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;
            let mut table_writer = table_path_to_writer(table);
            let mut seen = HashMap::<u64, String, BuildNoHashHasher<u64>>::default();

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, canonicalized| {
                    // runs in worker
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    *canonicalized = circkit::canonicalize(&normalized);
                },
                |record, canonicalized| {
                    // runs in main thread

                    let canonicalized_hash = xxhash_rust::xxh3::xxh3_64(canonicalized);

                    if !seen.contains_key(&canonicalized_hash) {
                        seen.insert(canonicalized_hash, record.id().unwrap().to_owned());

                        writer.write_all(b">").unwrap();
                        writer.write_all(record.head()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        match canonicalize {
                            true => {
                                writer.write_all(canonicalized).unwrap();
                            }
                            false => {
                                writer.write_all(record.seq()).unwrap();
                            }
                        };
                        writer.write_all(b"\n").unwrap();
                    } else {
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer
                                .serialize(Row {
                                    id: seen.get(&canonicalized_hash).unwrap(),
                                    duplicate_id: record.id().unwrap(),
                                })
                                .expect("failed to serialize table row");
                        }
                    }

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
            if let Some(mut table_writer) = table_writer {
                table_writer.flush()?;
            }
        }
        _ => panic!("input command is not for uniq"),
    }
    Ok(())
}
