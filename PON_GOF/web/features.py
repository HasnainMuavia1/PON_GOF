from .all_function import *


def extract_feature(recor_df):

    fasta_records = recor_df

    dpc_features = []
    residue_features = []
    position_features = []
    redundancy_features = []
    housekeeping_features = []
    ip_features = []
    essential_genes_features = []
    haploinsufficient_features = []
    lethality_features = []
    compl_knockout_features = []
    pssm_features = []

    accessibility_features = []
    idr_features = []
    transmembrane_features = []
    gene_age_features = []
    graph_features = []
    duplicate_gene_feature = []
    pseudo_gene_features = []

    protDcal_features = []
    repeats_features = []

    for i, row in fasta_records.iterrows():
        print("Features Calculation for: ", i, row['refseq_ids'])
        refseq_ids = row['refseq_ids']
        var_id = row['variation_ids']
        seq = mane_seqs.get(refseq_ids)
        ensgene_id, gene_id, ensprot, nm_ids = data_dict.get(refseq_ids)
        uniprot_id = ref_dict_file.get(refseq_ids)
        ##DPC features
        dpc_df = cal_dipeptide(seq)
        dpc_features.append(dpc_df)
        #residue featue
        res_df = get_residue(seq, var_id)
        residue_features.append(res_df)
        #position featues
        first_pos_df = get_first_position(refseq_ids, var_id)
        position_features.append(first_pos_df)
        #evlutionary features
        pssm_feature_df = get_pssm_features(refseq_ids, var_id)
        pssm_features.append(pssm_feature_df)
        #IDR features
        idr_feature_df = get_IDR_region(uniprot_id, var_id, IDR_records)
        idr_features.append(idr_feature_df)
        #Transmembrane features
        transmembrane_feature_df = get_transmembran_region(uniprot_id, var_id, htp_data, refseq_ids)
        transmembrane_features.append(transmembrane_feature_df)
        #Gene age features
        gene_age_feature_df = get_gene_age(refseq_ids, gene_age_data)
        gene_age_features.append(gene_age_feature_df)

        data_record = data_dict.get(refseq_ids)
        ensprot_id = data_record[2]
        #Graph features
        graph_df = get_graph_features(ensprot_id, graph_data)
        graph_features.append(graph_df)
        #duplicate featues
        duplct_gene_df = get_duplicate_gene_feature(refseq_ids, duplct_gene_data)
        duplicate_gene_feature.append(duplct_gene_df)
        #pseudo gene feautures
        pseudo_gene_df = get_pseudogene_feature(gene_id, pseudogene_dict)
        pseudo_gene_features.append(pseudo_gene_df)
        #protDcal features
        thermo_dynaminc_df = get_thermo_dynaminc_feature(refseq_ids, var_id)
        protDcal_features.append(thermo_dynaminc_df)
        #repeats features
        repeats_feature_df = get_repeats_region(refseq_ids, var_id, repeats_data)
        repeats_features.append(repeats_feature_df)

        # Accessibility features
        if refseq_ids:
            accessibility_df = get_accessibility_feature(refseq_ids, var_id)
            accessibility_features.append(accessibility_df)
        # Gene based features
        if gene_id != 0:
            red_df = get_redundancy(gene_id, family_list, var_id)
            redundancy_features.append(red_df)
            hkf_df = get_housekeeping(gene_id, hk_records, var_id)
            housekeeping_features.append(hkf_df)
            ip_df = get_interitencePatterns(gene_id, InheritenceP_records)
            ip_features.append(ip_df)
            ess_df = get_essentialgenes(gene_id, ess_genes_records)
            essential_genes_features.append(ess_df)
            hsp_df = get_haploinsufficient_proteins(gene_id, HSP_records)
            haploinsufficient_features.append(hsp_df)
            lethal_df = get_lethality_features(gene_id, lethality_genes_list)
            lethality_features.append(lethal_df)
            ckf_df = get_comp_knockout_features(gene_id, ckout_genes_list)
            compl_knockout_features.append(ckf_df)

    # Convert lists of DataFrames into single DataFrames
    dpc_f_df = pd.concat(dpc_features, ignore_index=True)
    re_f_df = pd.concat(residue_features, ignore_index=True)
    pos_f_df = pd.concat(position_features, ignore_index=True)
    redundancy_f_df = pd.concat(redundancy_features, ignore_index=True)
    housekeeping_f_df = pd.concat(housekeeping_features, ignore_index=True)
    ip_f_df = pd.concat(ip_features, ignore_index=True)
    essential_genes_f_df = pd.concat(essential_genes_features, ignore_index=True)
    haploinsufficient_f_df = pd.concat(haploinsufficient_features, ignore_index=True)
    lethality_f_df = pd.concat(lethality_features, ignore_index=True)
    compl_knockout_f_df = pd.concat(compl_knockout_features, ignore_index=True)
    pssm_f_df = pd.concat(pssm_features, ignore_index=True)
    acesblty_f_df = pd.concat(accessibility_features, ignore_index=True)
    idr_f_df = pd.concat(idr_features, ignore_index=True)
    htp_f_df = pd.concat(transmembrane_features, ignore_index=True)
    age_f_df = pd.concat(gene_age_features, ignore_index=True)
    graph_f_df = pd.concat(graph_features, ignore_index=True)
    duplicate_f_df = pd.concat(duplicate_gene_feature, ignore_index=True)
    pseudogene_f_df = pd.concat(pseudo_gene_features, ignore_index=True)
    protDcal_f_df = pd.concat(protDcal_features, ignore_index=True)
    repeats_f_df = pd.concat(repeats_features, ignore_index=True)

    combined_features = pd.concat([dpc_f_df, re_f_df, pos_f_df,
                                   redundancy_f_df, housekeeping_f_df, ip_f_df,
                                   essential_genes_f_df, haploinsufficient_f_df,
                                   lethality_f_df, compl_knockout_f_df, pssm_f_df, acesblty_f_df, idr_f_df,
                                   htp_f_df, age_f_df, graph_f_df, duplicate_f_df, pseudogene_f_df, protDcal_f_df, repeats_f_df], axis=1)

    return combined_features
