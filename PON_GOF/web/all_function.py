from collections import defaultdict
import pandas as pd
import os, io
import re
import json
import numpy as np
from Bio import Entrez
from django.conf import settings

Entrez.email = 'saeed.ahmed@med.lu.se'

# Define paths


a_list = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
aa_list = [i + j for i in a_list for j in a_list]

# Read the family file and create a list of lists - This file is required for get_redundancy feature
# hs_family_file = wd+'data_files/hsapiens.Families.Strict.2R.txt'
hs_family_file = os.path.join(settings.DATA_FILES_DIR,'hsapiens.Families.Strict.2R.txt')
with open(hs_family_file, 'r') as f:
    family_list = [line.strip().split('\t') for line in f]

# Read the House Keeping file and create a list of lists
hk_records = pd.read_csv(os.path.join(settings.DATA_FILES_DIR,"Housekeeping_GenesHuman.csv"))

# Read the family file and create a list of lists
InheritenceP_records = pd.read_excel(os.path.join(settings.DATA_FILES_DIR,"Clinical_Genomic_Database.xls"))

# Read the family file and create a list of lists
ess_genes_records = pd.read_excel(os.path.join(settings.DATA_FILES_DIR,"DEG_Info_Essential_genes.xlsx"))

# Read the family file and create a list of lists
HSP_records = pd.read_excel(os.path.join(settings.DATA_FILES_DIR,"Haploinsufficient_proteins.xlsx"))

# Read the family file and create a list of lists
lethality_genes_list = pd.read_excel(os.path.join(settings.DATA_FILES_DIR,"lethality_genes.xlsx"))['Lethality_Genes'].tolist()

# Read the family file and create a list of lists
ckout_genes_list = pd.read_excel(os.path.join(settings.DATA_FILES_DIR,"Complete_knockout_genes.xlsx"))['knockout_genes'].tolist()

# Read the IDR file and create a list of lists
IDR_records = pd.read_csv(os.path.join(settings.DATA_FILES_DIR,"IDR_data.csv"))

# Read file for Transmembrane regions data - HTP Human Transmembrane Protein
htp_data = pd.read_csv(os.path.join(settings.DATA_FILES_DIR,"Htp_data.csv"))

# GENE Age data file
gene_age_data = pd.read_csv(os.path.join(settings.DATA_FILES_DIR,"gene_age_data.csv"))

graph_data = pd.read_csv(os.path.join(settings.DATA_FILES_DIR, 'Graph_features_measures.csv'))

# Duplicate Gene Data - Simplified File
duplct_gene_data = pd.read_csv(os.path.join(settings.DATA_FILES_DIR, 'duplicate_gene_data_simplified.csv'))

# Pseudogene data
pseudogene_file = pd.read_csv(os.path.join(settings.DATA_FILES_DIR, 'pseudogene_output.txt'), sep='\t')
pg_gene = pseudogene_file['gene_symbol']
pg_value = pseudogene_file['pseudo_gene_status']
pseudogene_dict = dict(zip(pg_gene, pg_value))

# Pfam dictionary
ensgene_pfam_file = pd.read_csv(os.path.join(settings.DATA_FILES_DIR, 'ensemblToPfam_dict.txt'), sep='\t')
ens_gene_ids = ensgene_pfam_file['Ensembl_Gene_ID']
pfam_ids = ensgene_pfam_file['Pfam']
pfam_data_dict = dict(zip(ens_gene_ids, pfam_ids))

# Dictionary for ISA-ECI features
data_file = os.path.join(settings.DATA_FILES_DIR, 'amino_acid_data.json')
with open(data_file, 'r') as file:
    dict_file = json.load(file)

# Dictionary for uniprot_refseq features
data_file = os.path.join(settings.DATA_FILES_DIR, 'uniport_to_refseq.json')
with open(data_file, 'r') as file:
    ref_dict_file = json.load(file)

# Dictionary for uniprot accession numbers for alphafold
gene_dict_file = os.path.join(settings.DATA_FILES_DIR, 'refseq_to_uniprot_alphafold_dict.txt')
refseq_uniprot_dict = {}
with open(gene_dict_file, 'r') as file:
    for line in file:
        line = line.strip()
        parts = line.split('\t')
        if len(parts) == 2:
            prot_id, gene_name = parts[0], parts[1]
            refseq_uniprot_dict[prot_id] = gene_name

print("Dictionary related to alphafold uniport ids are loaded successfully.")

# Dictionary for EnsGene to EnsProt
data_file = os.path.join(settings.DATA_FILES_DIR, 'dict_for_mane_to_other_systems.json')
with open(data_file, 'r') as file:
    data_dict = json.load(file)

print("Dictionary related to EnsGene and EnsProt are loaded successfully.")

# MANE sequences dictionary
seqs_file = os.path.join(settings.DATA_FILES_DIR, 'MANE_seqs_dict.json')
with open(seqs_file, 'r') as file:
    mane_seqs = json.load(file)

print("Dictionary related to MANE sequences loaded successfully.")


def check_aa(seq, aa):
    seq = "".join(seq.split())
    aaf = aa[0]  # from
    aai = int(aa[1:-1])  # index
    aat = aa[-1]  # to
    return seq, aaf, aat, aai


def msg_find(seq, aa):
    seq, aaf, aat, aai = check_aa(seq, aa)
    msg = ""
    if aaf not in a_list:
        msg = "aa error, origin of aa is invalid."
    if aat not in a_list:
        msg = "aa error, nutation of aa is invalid."
    if aai < 1 or aai > len(seq):
        msg = "aa error, index of aa is invalid."
    if seq[aai - 1] != aaf:
        msg = "aa error, seq[{}] = {}, but origin of aa = {}".format(aai, seq[aai - 1], aaf)


# Di-peptides features calculation
def cal_dipeptide(sequence):
    dipeptides = [a + b for a in a_list for b in a_list]
    dp_counts = {dp: 0 for dp in dipeptides}
    for i in range(len(sequence) - 1):
        dp = sequence[i:i + 2]
        if dp in dipeptides:
            dp_counts[dp] += 1
    total_count = sum(dp_counts.values())
    dp_comp = {dp: count / total_count for dp, count in dp_counts.items()}
    rounded_dp_comp = {dp: round(count, 4) for dp, count in dp_comp.items()}
    dpc_df = pd.DataFrame([rounded_dp_comp], columns=dipeptides)

    return dpc_df


def get_residue(seq, aa):
    seq, aaf, aat, aai = check_aa(seq, aa)
    res = aai
    res_df = pd.DataFrame({'Residue': [res]})
    return res_df


def get_first_position(seq, aa):
    seq, aaf, aat, aai = check_aa(seq, aa)
    pos = 0
    if aai == 1:
        pos = 1
    pos_df = pd.DataFrame({'First_pos': [pos]})
    return pos_df


def get_redundancy(gene_name, family_list, variation):
    family_id = 0
    for family in family_list:
        if gene_name in family:
            family_id = 1
            break
    family_id_df = pd.DataFrame({'Redundancy': [family_id]})
    return family_id_df


def get_housekeeping(geneID, hk_records, variation):
    hk_gene_ID = 0
    for i, row in hk_records.iterrows():
        gene_record = row['gene_name']
        if geneID == gene_record:
            hk_gene_ID = 1
            break
    hk_df = pd.DataFrame({'Housekeeping': [hk_gene_ID]})
    return hk_df


def get_interitencePatterns(geneID, IP_records):
    inheritance_counts = {'AD_iN': 0, 'AR_iN': 0, 'AD/AR': 0, 'XL': 0}
    result_list = []
    for i, row in IP_records.iterrows():
        if geneID in row['GENE']:
            inherit_record = row['INHERITANCE'].strip()
            if inherit_record in inheritance_counts:
                inheritance_counts[inherit_record] += 1
    result_list.append(inheritance_counts.copy())
    df = pd.DataFrame(result_list)
    return df


def get_essentialgenes(geneID, eg_records):
    e_gene = 0
    for i, row in eg_records.iterrows():
        es_gene = row['GENE']
        if geneID == es_gene:
            e_gene = 1
            break
    e_gene_df = pd.DataFrame({'Essential_Gene': [e_gene]})
    return e_gene_df


def get_haploinsufficient_proteins(geneID, HSP_records):
    hsp_value = 0
    for i, row in HSP_records.iterrows():
        hsp_gene = row['Gene Symbol']
        if geneID == hsp_gene:
            hsp_value = 1
            break
    hsp_gene_df = pd.DataFrame({'HSP_value': [hsp_value]})
    return hsp_gene_df


def get_lethality_features(geneID, lethal_list):
    lethality_gene = 0
    if geneID in lethal_list:
        lethality_gene = 1
    lethality_gene_df = pd.DataFrame({'Lethality': [lethality_gene]})
    return lethality_gene_df


def get_comp_knockout_features(geneID, ckout_list):
    knockout_gene = 0
    if geneID in ckout_list:
        knockout_gene = 1
    knockout_gene_df = pd.DataFrame({'Knockout_Gene': [knockout_gene]})
    return knockout_gene_df


# Evolutionary information - PSSM input file reading functions
def updated_read_pssm_matrix(input_matrix):
    """ Reading the PSSM input file into a numpy array"""
    PSSM = []
    p = re.compile(r'-*[0-9]+')
    stream = open(input_matrix)
    for line, string in enumerate(stream.readlines()):
        if line > 2:
            str_vec = []
            overall_vec = string.split()
            if len(overall_vec) == 0:
                break
            str_vec.extend(overall_vec[1])

            if len(overall_vec) < 44:
                for cur_str in overall_vec[2:]:
                    str_vec.extend(p.findall(cur_str))
                    if len(str_vec) >= 21:
                        if len(str_vec) > 21:
                            raise ValueError("Wrong PSSM format")
                        break
                print("Done")
            else:
                str_vec = string.split()[1:44]
            if len(str_vec) == 0:
                break
            PSSM.append(str_vec)
    PSSM = np.array(PSSM)
    return PSSM

def get_pssm_features(prot_id, var_id):
    wd_feat = os.path.join(settings.DATA_FILES_DIR, '09_pssm_files')
    match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
    if match:
        org_aa, loc_aa, var_aa = match.group(1), int(match.group(2)) - 1, match.group(3)  # -1 is because of 0-index

    pssm_file_path = os.path.join(wd_feat, f"{prot_id}.pssm")
    input_matrix = updated_read_pssm_matrix(pssm_file_path)
    column_names = ['aa', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                    'V', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                    'V', 'Last1', 'Last2']
    input_df = pd.DataFrame(input_matrix, columns=column_names)
    row_record = input_df.iloc[loc_aa]
    var_values = row_record[var_aa]
    value1, value2 = var_values
    value3, value4 = float(row_record['Last1']), float(row_record['Last2'])
    count_file = os.path.join(settings.DATA_FILES_DIR, 'all_alignments_count_updated', f"{prot_id}_counts.csv")
    count_read = pd.read_csv(count_file)
    count_df = pd.DataFrame(count_read)
    loc_aa = loc_aa + 1
    count_record = count_df[count_df['pos'] == loc_aa]
    count_value = count_record['ratio'].values
    count_value = count_value[0]

    # Creating a DataFrame with the PSSM features
    pssm_df = pd.DataFrame({'pssm1': [int(value1)],
                            'pssm2': [int(value2)],
                            'pssm3': [value3],
                            'pssm4': [value4],
                            'pssm5': [count_value]})
    return pssm_df


def get_accessibility_feature(uniprot_id, var_id):
    wd_feat = os.path.join(settings.DATA_FILES_DIR, '5_FreeSASA_accessibility_area_NP')
    match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
    if match:
        org_aa, loc_aa, var_aa = match.group(1), int(match.group(2)), match.group(3)

    key = f"{org_aa}_{loc_aa}"
    file_path = os.path.join(wd_feat, f"{uniprot_id}.txt")
    with open(file_path, 'r') as file:
        data = file.read()

    lines = data.split('\n')
    header = lines[0]
    data_lines = lines[1:]
    data = {}
    for line in data_lines:
        parts = line.strip().split(',')  # Assuming your data is comma-separated
        if len(parts) == 2:
            data[parts[0]] = parts[1]  # get the accessibility area value for each variation

    if key in data:
        value_str = data[key]  # get the corresponding new_area value for the variation
        value_list = eval(value_str)  # Convert the string representation of the list to an actual list

        headers = ['Accessibility']
        accessibility_df = pd.DataFrame([value_list], columns=headers)

        return accessibility_df


def get_IDR_region(prot_acc, var_id, dp_records):
    var_status = 0
    all_records = dp_records[dp_records['acc_ids'] == prot_acc]
    if all_records.empty:
        var_status = 0
    else:
        for ind, record in all_records.iterrows():
            # Assuming 'sequence' is a column in dp_records containing the protein sequence
            sequence = record['sequences']
            # Assuming 'start_position' and 'end_position' are columns in dp_records
            start_position = record['start_points']
            end_position = record['end_points']

            match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
            if match:
                org_AA = match.group(1)
                loc_AA = int(match.group(2))
                var_AA = match.group(3)
            if loc_AA >= start_position and loc_AA <= end_position:
                target_pos = loc_AA - start_position
                amino_acid = sequence[target_pos]
                if amino_acid == org_AA:
                    var_status = 1
                else:
                    var_status = 0
            else:
                var_status = 0

    headers = ['IDR_Feat']
    idr_df = pd.DataFrame([var_status], columns=headers)
    return idr_df


def check_variation_status(prot_sequence, var_id, prot_id):
    status = False
    if isinstance(var_id, str):
        match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
        if match:
            org_AA = match.group(1)
            loc_AA = int(match.group(2))
            var_AA = match.group(3)
            if loc_AA <= len(prot_sequence):
                amino_acid = prot_sequence[loc_AA - 1]
                if amino_acid == org_AA:
                    status = True
                else:
                    print(
                        f" Variation: {var_id}, '{org_AA}' not found at  Location {loc_AA} instead {amino_acid} is at Location {loc_AA}: ")
    else:
        print(f"Skipping invalid var_id: {prot_id} for {var_id}")
    return status


def get_transmembran_region(prot_acc, var_id, dp_records, prot_id):
    var_status = 0
    all_records = dp_records[dp_records['Accession'] == prot_acc]
    if all_records.empty:
        var_status = 0
    else:
        for ind, record in all_records.iterrows():
            # Assuming 'sequence' is a column in dp_records containing the protein sequence
            sequence = record['Sequence']
            var_status = check_variation_status(sequence, var_id, prot_id)
            if var_status == True:
                var_status = 1
            else:
                var_status = 0
    headers = ['HTP_Feat']
    htp_df = pd.DataFrame([var_status], columns=headers)
    return htp_df


def get_gene_age(refseq_id, gene_age_file):
    gene_age_value = 0
    age_data = gene_age_file[gene_age_file['refseq_ids'] == refseq_id]
    if not age_data.empty:
        gene_age_value = age_data['gene_age'].iloc[0]
    headers = ['Gene_age']
    gene_age_df = pd.DataFrame([gene_age_value], columns=headers)
    return gene_age_df


def get_graph_features(ens_gene_id, graph_data):
    graph_features = ['0', '0', '0', '0', '0', '0', '0', '0', '0']
    records_data = graph_data[graph_data['Proteins'] == ens_gene_id]
    if not records_data.empty:
        graph_features = records_data.iloc[:, 1:].values.flatten().tolist()
    else:
        graph_features

    headers = ['Degree', 'Closeness', 'Betweenness', 'Eigenvector', 'Harmonic', 'Hub_Score', 'Authority_Score',
               'Page_Rank', 'Power_Centrality']
    graph_features_df = pd.DataFrame([graph_features], columns=headers)
    return graph_features_df


def get_duplicate_gene_feature(refseq_id, data):
    data_list = data.values.flatten()
    data_record = data_dict.get(refseq_id)
    mRNA_id = data_record[3]
    mRNA_id = mRNA_id.split('.')[0]

    if mRNA_id in data_list:
        duplct_gene_value = 1
    else:
        duplct_gene_value = 0
    headers = ['Duplicate_Gene']
    duplct_gene_df = pd.DataFrame([duplct_gene_value], columns=headers)
    return duplct_gene_df


def get_pseudogene_feature(gene_id, pseudogene_dict):
    pg_gene_feature = 0
    gene_value_in_dict = pseudogene_dict.get(gene_id)
    if gene_value_in_dict:
        pg_gene_feature = gene_value_in_dict
    else:
        pg_gene_feature
    header = ['Pseudogene']
    pg_gene_df = pd.DataFrame([pg_gene_feature], columns=header)
    return pg_gene_df


def get_pfam_feature(gene_id, pseudogene_dict):
    pg_gene_feature = 0
    gene_value_in_dict = pseudogene_dict.get(gene_id)
    if gene_value_in_dict:
        pg_gene_feature = 1
    else:
        pg_gene_feature = 0
    header = ['pfam']
    pg_gene_df = pd.DataFrame([pg_gene_feature], columns=header)
    return pg_gene_df


def get_thermo_dynaminc_feature(refseq_id, var_id):
    # wd_feat = "/data_files/ProtDCal_MANE/"
    match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
    if match:
        org_aa, loc_aa, var_aa = match.group(1), int(match.group(2)), match.group(3)

    key = f"{org_aa}_{loc_aa}"
    file_path = os.path.join(settings.DATA_FILES_DIR, f"ProtDCal_MANE/{refseq_id}.txt")
    with open(file_path, 'r') as file:
        data = file.read()

    lines = data.split('\n')
    header = lines[0]
    data_lines = lines[1:]
    data_dict = {}
    for line in data_lines:
        parts = line.strip().split()
        if len(parts) == 4:
            data_dict[parts[0]] = parts[1:]

    if key in data_dict:
        value_list = data_dict[key]
        headers = ['Gw_U', 'Gs_U', 'W_U']
        thermo_dyanamic_fea_df = pd.DataFrame([value_list], columns=headers)
        return thermo_dyanamic_fea_df
    else:
        print("Key not found in data")
        return None


# Read file for Repeats data - Read only the specified columns from the repeats data file
repeats_data = pd.read_csv(os.path.join(settings.DATA_FILES_DIR,"repeat.csv"), usecols=['protein_id', 'from', 'to'])


def get_repeats_region(prot_id, var_id, repeats_records):
    match = re.match(r'([A-Z])(\d+)([A-Z])', var_id)
    if match:
        var_loc = int(match.group(2))

    all_records = repeats_records[repeats_records['protein_id'] == prot_id]
    var_status = 0  # Default status is 0 (False)
    for ind, record in all_records.iterrows():
        start = record['from']
        end = record['to']
        if start <= var_loc <= end:
            var_status = 1
            break

    headers = ['Repeats']
    repeats_df = pd.DataFrame([var_status], columns=headers)
    return repeats_df