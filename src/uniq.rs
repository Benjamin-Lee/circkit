use crate::{
    commands::Command,
    utils::{
        canonicalization_input, normalized_sequence, output_to_writer, process_fasta,
        table_path_to_writer,
    },
};
use anyhow::Context;
use nohash_hasher::BuildNoHashHasher;
use seq_io::fasta::Record;
use std::collections::{hash_map::Entry, HashMap};
use std::io::Write;

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
            table_format,
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
            let mut table_writer = table_path_to_writer(table, *table_format)?;
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
                            entry.insert(if table_writer.is_some() {
                                record
                                    .id()
                                    .context(
                                        "FASTA identifiers must be UTF-8 when writing metadata",
                                    )?
                                    .to_owned()
                            } else {
                                String::new()
                            });

                            writer.write_all(b">")?;
                            writer.write_all(record.head())?;
                            writer.write_all(b"\n")?;
                            match canonicalize {
                                true => {
                                    writer.write_all(&data.0)?;
                                }
                                false => {
                                    writer.write_all(record.seq())?;
                                }
                            };
                            writer.write_all(b"\n")?;
                        }
                        Entry::Occupied(entry) => {
                            if let Some(ref mut table_writer) = table_writer {
                                table_writer.serialize(Row {
                                    id: entry.get(),
                                    duplicate_id: record.id().context(
                                        "FASTA identifiers must be UTF-8 when writing metadata",
                                    )?,
                                })?;
                            }
                        }
                    }

                    Ok(())
                },
            )?;
            writer.finish()?;
            if let Some(table_writer) = table_writer {
                table_writer.finish()?;
            }
        }
        _ => panic!("input command is not for uniq"),
    }
    Ok(())
}
