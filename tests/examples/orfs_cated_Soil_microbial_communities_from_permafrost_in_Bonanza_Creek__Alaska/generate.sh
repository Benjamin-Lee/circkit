xz -d -k -c ../nim_cated_Soil_microbial_communities_from_permafrost_in_Bonanza_Creek__Alaska/in.fasta.xz > in.fasta
orfipy in.fasta --min 75  --dna out.fasta --strand f --include-stop
rm in.fasta
mv orfipy_in.fasta_out/out.fasta out.fasta
rm -r orfipy_in.fasta_out 