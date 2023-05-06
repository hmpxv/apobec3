

squirrel Clade_IIb_2022-10-06_with_gisaid.fasta


iqtree -asr -s Clade_IIb_2022-10-06_with_gisaid.aln.fasta -o "KJ642615|W-Nigeria|Nigeria|1978" -czb -m JC -T 1 

jclusterfunk prune -i  Clade_IIb_2022-10-06_with_gisaid.aln.fasta.treefile -o Clade_IIb_2022-10-06_with_gisaid.pruned.tree -t "KJ642615|W-Nigeria|Nigeria|1978"

