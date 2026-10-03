use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn tuning_preserves_fasta_and_metadata_in_both_execution_modes() {
    let directory = assert_fs::TempDir::new().unwrap();
    let input = directory.path().join("input.fasta");
    std::fs::write(
        &input,
        b">one description\nATGAAATAAATGAAATAA\n>two\natgcNN---\n>three\nCACACACACACA\n",
    )
    .unwrap();
    for operation in ["canonicalize", "uniq", "monomerize", "orfs"] {
        let reference = directory.path().join("reference.fasta");
        let reference_table = directory.path().join("reference.tsv");
        let mut base = Command::cargo_bin("circkit").unwrap();
        base.args([operation, "--threads", "1", "--execution", "serial"])
            .arg(&input)
            .arg("-o")
            .arg(&reference);
        if operation != "canonicalize" {
            base.arg("--table").arg(&reference_table);
        }
        if operation == "monomerize" {
            base.args(["--keep-all", "--max-mismatch", "2"]);
        }
        base.assert().success();
        for mode in ["serial", "pipeline"] {
            let output = directory.path().join("actual.fasta");
            let table = directory.path().join("actual.tsv");
            let mut tuned = Command::cargo_bin("circkit").unwrap();
            tuned
                .args([
                    operation,
                    "--threads",
                    "1",
                    "--execution",
                    mode,
                    "--queue-depth",
                    "1",
                ])
                .arg(&input)
                .arg("-o")
                .arg(&output);
            if operation != "canonicalize" {
                tuned.arg("--table").arg(&table);
            }
            if operation == "canonicalize" || operation == "uniq" {
                tuned.args(["--rotation-cutoff", "0"]);
            }
            if operation == "monomerize" {
                tuned.args([
                    "--keep-all",
                    "--max-mismatch",
                    "2",
                    "--mismatch-chunk-size",
                    "1",
                ]);
            }
            tuned.assert().success();
            assert_eq!(
                std::fs::read(&reference).unwrap(),
                std::fs::read(output).unwrap(),
                "{operation} {mode}"
            );
            if operation != "canonicalize" {
                assert_eq!(
                    std::fs::read(&reference_table).unwrap(),
                    std::fs::read(table).unwrap(),
                    "{operation} {mode} metadata"
                );
            }
        }
    }
}

#[test]
fn invalid_tuning_is_rejected_and_legacy_queue_flag_works() {
    for (operation, flag) in [
        ("orfs", "--threads"),
        ("canonicalize", "--queue-depth"),
        ("monomerize", "--mismatch-chunk-size"),
    ] {
        Command::cargo_bin("circkit")
            .unwrap()
            .args([operation, flag, "0"])
            .assert()
            .failure();
    }
    Command::cargo_bin("circkit")
        .unwrap()
        .args(["canonicalize", "--threads", "2", "--execution", "serial"])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Serial execution requires --threads 1",
        ));
    Command::cargo_bin("circkit")
        .unwrap()
        .args([
            "monomerize",
            "--threads",
            "1",
            "--batch-size",
            "1",
            "--keep-all",
        ])
        .write_stdin(">one\nACGT\n")
        .assert()
        .success()
        .stdout(">one\nACGT\n");
}
