import pandas as pd
import ast

# Input file paths
file1 = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Small_Experiment/Intepro_Analysis/Up/Up_IPR_Enrichment.FUNC-E.efeatures.tsv"
file2 = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Small_Experiment/Intepro_Analysis/Up/Up_IPR_Enrichment.FUNC-E.enriched_terms.tsv"

# Load gene annotations
df_genes = pd.read_csv(file1, sep="\t")
df_genes["Term"] = df_genes["Term"].apply(ast.literal_eval)

# Load enriched IPR terms
df_enriched = pd.read_csv(file2, sep="\t")

# Prepare output list
output_rows = []

# For each enriched IPR term, find matching genes in df_genes
for _, enrich_row in df_enriched.iterrows():
    ipr = enrich_row["Term"]
    name = enrich_row["Name"]
    
    for _, gene_row in df_genes.iterrows():
        gene = gene_row["Feature"]
        gene_iprs = gene_row["Term"]
        
        if ipr in gene_iprs:
            output_rows.append({
                "IPR Term": ipr,
                "Name": name,
                "Gene ID": gene
            })

# Convert to DataFrame and save
df_output = pd.DataFrame(output_rows)
output_path = "/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Small_Experiment/Intepro_Analysis/Up/Up_Enriched_Gene_To_IPR.tsv"
df_output.to_csv(output_path, sep="\t", index=False)

print(f"✅ Output written to: {output_path}")

