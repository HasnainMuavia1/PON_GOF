import os
import pandas as pd
import numpy as np
import joblib
from flask import Flask, request, redirect, flash, render_template
from features import extract_feature
from tabulate import tabulate

app = Flask(__name__)
app.secret_key = "supersecretkey"

# Define paths

genes_dict_path = '/data_files/Updated_Ens_genes_dict.txt'
protein_gene_info = {}
with open(genes_dict_path, 'r') as file:
    for line in file:
        line = line.strip()
        parts = line.split('\t')
        if len(parts) == 2:
            prot_id, gene_name = parts[0], parts[1]
            protein_gene_info[prot_id] = gene_name
print("Dictionary loaded successfully.")

go_file = pd.read_csv('/data_files/gene_GO_sequence.csv')
go_ratio_dict = dict(zip(go_file['ensembl_gene_ids'], go_file['GO_ratio']))

models_dir = os.path.join("/models")
select_indices_file_path = os.path.join('selected_feat_indecs_50.csv')


@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file:
            file_path = os.path.join('input_file.csv')
            file.save(file_path)
            predictions, headers = getPrediction(file_path)
            formatted_table = tabulate(predictions, headers=headers, floatfmt=".3f")
            return render_template('index.html', predictions=formatted_table)
    return render_template('index.html', predictions=None)


@app.route('/clear', methods=['GET'])
def clear():
    return redirect('/')


def getPrediction(file_name):
    record_file = pd.read_csv(file_name)
    test_feature = extract_feature(record_file)
    select_indices_file = pd.read_csv(select_indices_file_path)
    selected_features = select_indices_file['Feature'].tolist()
    test_data = test_feature[selected_features]

    fs_ids = record_file['refseq_ids']
    refseq_ids = []
    go_ratio_value = []

    for id in fs_ids:
        ensembl_id = protein_gene_info.get(id)
        go_value = go_ratio_dict.get(ensembl_id, 0)
        refseq_ids.append(id)
        go_ratio_value.append(go_value)

    go_data = pd.DataFrame({'refseq_ids': refseq_ids, 'GO_ratio': go_ratio_value})
    test_data['GO_value'] = pd.to_numeric(go_data['GO_ratio'])
    test_data_wGO = test_data.values

    n_iterations = 200
    pred_matrix = np.empty((len(test_data_wGO), 0))

    for i in range(n_iterations):
        model_filename = os.path.join(models_dir, f'lgbm_model_{i + 1}.joblib')
        lgbm = joblib.load(model_filename)
        pred_probs = lgbm.predict_proba(test_data_wGO)[:, 1]
        pred_matrix = np.column_stack([pred_matrix, pred_probs])

    mean_probs = np.mean(pred_matrix, axis=1)
    pred_labels = np.where(mean_probs >= 0.5, 'P', 'N')

    results_table = pd.DataFrame({
        "refseq_ids": record_file["refseq_ids"],
        "variation_ids": record_file["variation_ids"],
        "meanProb": np.round(mean_probs, 3),
        "pred_class": pred_labels,
    })

    table_headers = ["refseq_ids", "variation_ids", "meanProb", "pred_class"]
    table_data = results_table.values.tolist()

    # Print results in table format for debugging
    print(tabulate(table_data, headers=table_headers, floatfmt=".3f"))

    return table_data, table_headers


if __name__ == "__main__":
    app.run(debug=True)
