import csv
from Bio import SeqIO

## input files required (called whatever)
reconstructed_branch_snps = "Clade_IIb_2022-08-22_with_gisaid.aln.pruned.tree.branch_snps.reconstruction.csv"
my_ref_file = '/Users/s1680070/repositories/alignHPXV/squirrel/data/NC_063383.fasta'
alignment_file = "Clade_IIb_2022-08-22_with_gisaid.aln.fasta"

#files created, name them whatever you like
apobec_mutation_file = "partition_analysis/apobec3_mutations.tsv"
apobec_partition_file = "partition_analysis/apobec3_only.fasta"
non_apobec_partition_file = "partition_analysis/non_apobec3_only.fasta"
unmasked_alignment_file = "partition_analysis/Clade_IIb_2022-08-22.unmasked.fasta"
apobec_index_record_file = "./updated_2022-09-02/partition_analysis/apobec3_constant_target_modified.csv"

apobec_variable = collections.defaultdict(list)
non_apobec_variable =  collections.defaultdict(list)

with open(reconstructed_branch_snps,"r") as f:
    reader = csv.DictReader(f)
    for row in reader:
#         print(row)
        
        if row["dimer"] in ["GA","TC"]:
            
            apobec_variable[int(row['site'])].append(f"{row['parent']}|{row['child']}")
        else:
            if row["site"] != "3450":
                non_apobec_variable[int(row['site'])].append(row)
            
print(len(apobec_variable), len(non_apobec_variable))

ref = str(SeqIO.read(my_ref_file, 'fasta').seq)
pos = np.arange(len(ref))

apo_keep_0index = set()
for i in pos:
    if ref[i:i+2]=="GA":
        if i+1 not in non_apobec_variable:
            apo_keep_0index.add(i)
    elif ref[i:i+2]=="TC":
        if i+2 not in non_apobec_variable:
            apo_keep_0index.add(i+1)
        
for i in apobec_variable:
    apo_keep_0index.add(i-1)

ftrait = open(apobec_mutation_file,"w")

ftrait.write("\t")
muts = []
for i in sorted(apobec_variable):
    muts.append(str(i))
muts = "\t".join(muts)    
ftrait.write(f"{muts}\n")

with open(apobec_partition_file,"w") as fw:
    with open(non_apobec_partition_file,"w") as fw2:
        with open(unmasked_alignment_file,"w") as fw3:
            for record in SeqIO.parse(alignment_file,"fasta"):
                
                new_id = record.id
                if new_id in names:
                    new_id = names[new_id]
                    
                if not record.id == "KJ642615|W-Nigeria|Nigeria|1978":
                    
                    ftrait.write(f"{new_id}\t")
                    
                    
                    muts = []
                    for i in sorted(apobec_variable):
                        if record.seq[i-1] in ["A","T"]:
                            muts.append(record.seq[i-1])
                        else:
                            muts.append(record.seq[i-1])
#                     print(muts)
                    muts = "\t".join(muts)    
                    ftrait.write(f"{muts}\n")
                    
                    apo_seq = list(str(record.seq))
                    non_apo_seq = list(str(record.seq))
                    
                    if record.id == "MT250197|Singapore_2019|Singapore_ex_Nigeria|2019-05-08":
                        apo_seq[3449] = "N"
                    apo_masked = 0
                    non_apo_masked = 0
                    for i in range(len(apo_seq)):
                        if i not in apo_keep_0index:
                            non_apo_masked +=1
                            non_apo_seq[i] = "N"
                        else:
                            apo_masked +=1
                            apo_seq[i] = "N"
#                     print(apo_masked,non_apo_masked)
                    apo_seq = "".join(apo_seq)
                    fw2.write(f">{new_id}\n{apo_seq}\n")
                    
                    non_apo_seq = "".join(non_apo_seq)
                    fw.write(f">{new_id}\n{non_apo_seq}\n")
                
                fw3.write(f">{new_id}\n{record.seq}\n")
            
            
with open(apobec_index_record_file,"w") as fw:
    fw.write("index\n")
    for i in apo_keep_0index:
        fw.write(f"{i}\n")
ftrait.close()           
        