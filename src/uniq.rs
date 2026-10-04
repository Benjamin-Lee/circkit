use crate::{
    commands::Command,
    utils::{
        canonicalization_input, normalized_sequence, output_to_writer, process_fasta,
        table_path_to_writer,
    },
};
use nohash_hasher::BuildNoHashHasher;
use seq_io::fasta::Record;
use std::collections::{hash_map::Entry, HashMap};

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
            processing,
            rotation,
        } => {
            let (reader, settings) = canonicalization_input(input, output, processing, *threads)?;
            let reuse_buffers = settings.serial;
            let rotation = circkit::canonicalize::RotationOptions {
                duval_max_len: rotation.rotation_cutoff,
            };
            let mut writer = output_to_writer(output)?;
            let mut table_writer = table_path_to_writer(table);
            let mut seen = HashMap::<u64, String, BuildNoHashHasher<u64>>::default();

            process_fasta(
                reader,
                settings,
                |record, data: &mut (Vec<u8>, Vec<u8>)| {
                    // runs in worker
                    let normalized = normalized_sequence(record.seq());
                    if reuse_buffers {
                        rotation.canonicalize_into(&normalized, &mut data.0, &mut data.1);
                    } else {
                        data.0 = rotation.canonicalize(&normalized);
                        data.1 = Vec::new();
                    }
                },
                |record, data| {
                    // runs in main thread

                    let canonicalized_hash = xxhash_rust::xxh3::xxh3_64(&data.0);

                    match seen.entry(canonicalized_hash) {
                        Entry::Vacant(entry) => {
                            entry.insert(record.id().unwrap().to_owned());

                            writer.write_all(b">").unwrap();
                            writer.write_all(record.head()).unwrap();
                            writer.write_all(b"\n").unwrap();
                            match canonicalize {
                                true => {
                                    writer.write_all(&data.0).unwrap();
                                }
                                false => {
                                    writer.write_all(record.seq()).unwrap();
                                }
                            };
                            writer.write_all(b"\n").unwrap();
                        }
                        Entry::Occupied(entry) => {
                            if let Some(ref mut table_writer) = table_writer {
                                table_writer
                                    .serialize(Row {
                                        id: entry.get(),
                                        duplicate_id: record.id().unwrap(),
                                    })
                                    .expect("failed to serialize table row");
                            }
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
