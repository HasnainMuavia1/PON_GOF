import os
import json

# Define the working directory
wd = "/media/saeedpsb/data/PON-GoF/Update_1/"

# Define the path to the input text file
input_file_path = os.path.join(wd, 'MANE_Select', 'MANE.GRCh38.v1.3.summary.txt')

# Initialize an empty dictionary
genes_dict = {}

try:
    # Open the input text file and process each line
    with open(input_file_path, 'r') as file:
        # Skip the header line
        next(file)

        # Process each line in the file
        for line in file:
            # Strip leading/trailing whitespace characters
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Split the line by tab character
            parts = line.split('\t')

            # Ensure the line has the correct number of columns
            if len(parts) < 14:
                print(f"Ignoring line with incorrect format: {line}")
                continue

            # Extract the RefSeq_prot and symbol fields
            refseq_prot = parts[6]
            symbol = parts[3]

            # Add the RefSeq_prot and symbol to the dictionary
            genes_dict[refseq_prot] = symbol
except FileNotFoundError:
    print(f"File not found: {input_file_path}")
except Exception as e:
    print(f"An error occurred: {e}")

# Define the path to the output dictionary file
output_file_path = os.path.join(wd, 'data_files', 'refseq2gensymob.txt')

try:
    # Write the dictionary to the output text file
    with open(output_file_path, 'w') as file:
        for key, value in genes_dict.items():
            file.write(f"{key}\t{value}\n")
except Exception as e:
    print(f"An error occurred while writing to the file: {e}")

# Print the dictionary to check the content (optional)
print(genes_dict)
