import pandas as pd

df = pd.read_csv("refseq_to_uniprot_update.csv")
uniport_to_ref_dict = dict(zip(df['UniProt ID'], df['RefSeq ID']))
# Save DataFrame to CSV file
df.to_csv('uniport_to_refseq.csv', index=False)

idr_df = pd.read_csv('IDR_data.csv')

import pandas as pd

# Step 1: Read the UniProt to RefSeq mapping CSV into a dictionary
refseq_to_uniprot_file = "uniport_to_refseq.csv"
df_mapping = pd.read_csv(refseq_to_uniprot_file)

# Create dictionary mapping UniProt IDs to RefSeq IDs
uniport_to_ref_dict = dict(zip(df_mapping['UniProt ID'], df_mapping['RefSeq ID']))

# Step 2: Read the IDR data CSV
idr_data_file = "IDR_data.csv"
idr_df = pd.read_csv(idr_data_file)

# Step 3: Extract RefSeq IDs based on UniProt IDs
idr_df['RefSeq ID'] = idr_df['acc_ids'].apply(lambda x: uniport_to_ref_dict.get(x, 'Not Found'))

# Step 4: Reorder columns to have RefSeq ID as the first column
idr_df = idr_df[['RefSeq ID'] + [col for col in idr_df.columns if col != 'RefSeq ID']]

# Step 5: Write the updated DataFrame back to CSV
idr_df.to_csv("IDR_data_with_RefSeq.csv", index=False)

print(f"Updated CSV file 'IDR_data_with_RefSeq.csv' has been created successfully.")

