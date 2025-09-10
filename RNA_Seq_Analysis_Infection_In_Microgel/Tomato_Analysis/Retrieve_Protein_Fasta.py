from Bio import SeqIO

# File paths
fasta_file = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Tomato_Genom_ITAG4/ITAG4.0_proteins.fasta"
ids_file = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis/Tip_VS_Root/All_Up.txt"
output_file = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis/Tip_VS_Root/All_Up_Proteins.fasta"

# Load target IDs into a set
with open(ids_file) as f:
    target_ids = set(line.strip() for line in f if line.strip())

# Extract matching sequences
count = 0
with open(output_file, "w") as out_f:
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Match full header or prefix (before .1/.2 etc.)
        if record.id in target_ids or record.id.split('.')[0] in target_ids:
            SeqIO.write(record, out_f, "fasta")
            count += 1

print(f"✅ Extracted {count} protein sequences into '{output_file}'.")

