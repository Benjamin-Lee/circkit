hyperfine --warmup 3 \
    --prepare 'rm -f tmp.fasta' \
    './target/release/circkit orfs data/img.monomers95.fasta --start-codons ATG,CTG,TTG --max-wraps 0 --include-stop --strand both -o tmp.fasta' \
    --prepare 'rm -f orfipy_img.monomers95.fasta_out/tmp.fasta' \
    'orfipy --min 75 data/img.monomers95.fasta --dna tmp.fasta --strand b --include-stop'