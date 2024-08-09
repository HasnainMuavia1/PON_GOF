'''
This code make a dictionary for refseq ids to accession numbers
used by Alpha fold
'''
import pandas as pd


csv_file_path = pd.read_csv("/data_files/ENSGene_to_ENSProt.csv")
ensgene = csv_file_path['ENSGene']
ensprot = csv_file_path['ENSProt']
ensgene_to_ensprot_dict = dict(zip(ensgene, ensprot))

# Write the dictionary to a file
gene_dict_file = "/data_files/ENSGene_to_ENSProt_dict.txt"
with open(gene_dict_file, "w") as f:
    for key, value in ensgene_to_ensprot_dict.items():
        f.write(f"{key}\t{value}\n")