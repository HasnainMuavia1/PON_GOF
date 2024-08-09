import pandas as pd
import json

# Read CSV into DataFrame
df = pd.read_csv("refseq_to_uniprot_update.csv")

# Create dictionary from DataFrame columns
uniport_to_ref_dict = dict(zip(df['UniProt ID'], df['RefSeq ID']))

# Save dictionary as JSON file
with open('uniport_to_refseq.json', 'w') as json_file:
    json.dump(uniport_to_ref_dict, json_file)

print("JSON file saved successfully.")

