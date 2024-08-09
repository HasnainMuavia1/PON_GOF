import pandas as pd

csv_file_path = pd.read_csv("/data_files/duplicate_gene_data.csv")
cluster_heads = csv_file_path['ClusterHead'].tolist()
members = csv_file_path['Member'].tolist()

# Combine the lists and remove duplicates while maintaining order
merged_data = list(dict.fromkeys(cluster_heads + members))

# Write the list to a file
duplicate_file = "/data_files/duplicate_gene_data_simplified.csv"
with open(duplicate_file, "w") as f:
    for item in merged_data:
        f.write("%s\n" % item)
