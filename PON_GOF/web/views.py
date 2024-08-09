import io

from django.shortcuts import render, HttpResponse, redirect
from .models import data,Counter
import csv
import pandas as pd
from django.core.mail import EmailMessage
import os
from django.db import transaction
from .features import extract_feature
import numpy as np
import joblib
import tabulate
from django.conf import settings


def index(request):
    return render(request, 'index.html')


# Create your views here.
def input(request):
    if request.method == 'POST':
        counter, created = Counter.objects.get_or_create(pk=1)  # Assuming a single counter instance
        counter.increment()
        print('hello')
        text = request.POST.get('te')
        file = request.FILES.get('csv_data')
        print(file)
        email = request.POST['e']
        print('hello2')

        if not text:
            # if the file is excel file then convert into csv
            if file.name.endswith('.xlsx') or file.name.endswith('.xls'):
                file = pd.read_excel(file)
                print(file.head())
                csv_file = file.to_csv(index=False)
            print('hello3')

            result = handle_uploaded_file2(file, email)
            if result == 'Error':
                return render(request, 'input_Prediction.html', {'msg': 'Wrong format or one of the keys is missing','counter':Counter.objects.first()})
            elif result =='Model error':
                return render(request, 'input_Prediction.html', {'msg': 'The values you provided are wrong ! unable to predict','counter':Counter.objects.first()})
            else:
                return render(request, 'input_Prediction.html', {'msg': 'Email sent successfully! Check your email','counter':Counter.objects.first()})
        elif not file:
            result = get_values_combine(text, email)
            if result == 'Error':
                return render(request, 'input_Prediction.html', {'msg': 'Failed to send email','counter':Counter.objects.first()})
            elif result =='Model error':
                return render(request, 'input_Prediction.html', {'msg': 'The values you provided are wrong ! unable to predict','counter':Counter.objects.first()})
            else:
                return render(request, 'input_Prediction.html', {'msg': 'Email sent successfully! Check your email','counter':Counter.objects.first()})

    return render(request, 'input_Prediction.html',{'counter':Counter.objects.first()})


def handle_uploaded_file2(file, email):
    """
    Process the uploaded CSV file and create a combined CSV file based on filtered data.
    """
    # Define required keys as a set
    required_keys = {'refseq_ids', 'variation_ids'}

    # Read the CSV file into a DataFrame
    df = pd.read_csv(file)

    # Get the DataFrame columns as a set
    actual_columns = set(df.columns)

    # Check if all required keys are in the actual columns
    if required_keys.issubset(actual_columns):
        # Group by 'refseq_ids'
        grouped = df.groupby('refseq_ids')
        df = pd.DataFrame()
        # List to store combined data and not-found data
        combined_data = []
        not_found_data = []

        # Process each group
        for var1, group in grouped:
            var2_list = group['variation_ids'].tolist()
            var2_list = [k.upper() for k in var2_list]
            print(var1)
            print(var2_list)

            # Query the database
            results = data.objects.filter(refseq_ids=var1, variation_ids__in=var2_list)

            # Convert the results to a list of dictionaries
            data_list = list(results.values('refseq_ids', 'variation_ids', 'meanProb', 'stdProb', 'pred_label'))
            combined_data.extend(data_list)

            # Track found variation_ids for the current group
            if data_list:
                found_variation_ids = {result['variation_ids'] for result in data_list}

                # Track not-found values for the current group
                not_found_variation_ids = [var2 for var2 in var2_list if var2 not in found_variation_ids]
            else:
                print("else condition owrkin")
                # If no data is found, consider all variations as not found
                not_found_variation_ids = var2_list

            # Append not-found data to the global list
            for var2 in not_found_variation_ids:
                not_found_data.append({
                    'refseq_ids': var1,
                    'variation_ids': var2,
                })
        print('before if condition')
        # Perform operations with not-found data only if not_found_data is not empty
        if not_found_data:
            not_found_df = pd.DataFrame(not_found_data)
            # file = not_found_df.to_csv(index=False)
            try:
                print('try condition')
                print("printing file head",not_found_df.head())
                df = get_prediction(not_found_df)
                if not df.empty:
                    model_predicted_csv = df.to_csv(index=False)
                    process_csv(model_predicted_csv)
                print('after process file ')
            except:
                print('model error')
        else:
            df = pd.DataFrame()  # Create an empty DataFrame if no data is found
        print('after else conditoin')
        # Convert the combined data list to a DataFrame
        combined_df = pd.DataFrame(combined_data)
        combined_df.drop_duplicates(subset=['refseq_ids', 'variation_ids'], inplace=True)

        # Combine the combined_df and not_found_df
        if not df.empty:
            combined_df = pd.concat([combined_df, df], ignore_index=True)

        # Create a CSV file from the combined DataFrame
        csv_buffer = io.StringIO()
        combined_df.to_csv(csv_buffer, index=False)
        csv_buffer.seek(0)  # Go to the start of the StringIO buffer

        # Send the email with CSV attachment
        email_message = EmailMessage(
            'Your files are here',
            'Hello, this is your file',
            to=[email]
        )
        email_message.attach('combined_data.csv', csv_buffer.getvalue(), 'text/csv')
        email_message.send()
        return 'success'
    else:
        # If the required keys are not present, return an error message
        return 'Error'


def get_values_combine(text, email):
    blocks = text.split('>')
    combined_data = []
    not_found_data = []
    df = pd.DataFrame()
    for block in blocks:
        block = block.strip()
        if not block:
            continue

        lines = block.split('\n')
        if lines:
            var1 = lines[0].strip()
            var2_list = [item.strip().upper() for item in lines[1:] if item.strip()]

            # Query the database
            results = data.objects.filter(refseq_ids=var1, variation_ids__in=var2_list)

            # Convert the results to a list of dictionaries
            data_list = list(results.values('refseq_ids', 'variation_ids', 'meanProb', 'stdProb', 'pred_label'))
            combined_data.extend(data_list)

            # Track found variation_ids for the current block
            if data_list:
                found_variation_ids = {result['variation_ids'] for result in data_list}
                print('found variation ids : ',found_variation_ids)

                # Track not found values for the current block
                not_found_variation_ids = [var2 for var2 in var2_list if var2 not in found_variation_ids]
                print('not found variation ids : ',not_found_variation_ids)

            else:
                # If no data is found, consider all variations as not found
                not_found_variation_ids = var2_list

            # Append not-found data to the global list
            for var2 in not_found_variation_ids:
                not_found_data.append({
                    'refseq_ids': var1,
                    'variation_ids': var2,
                })

    # Perform operations with not-found data only if not_found_data is not empty
    if not_found_data:
        not_found_df = pd.DataFrame(not_found_data)
        # file = not_found_df.to_csv(index=False)
        try:
            print('in value condition')
            df = get_prediction(not_found_df)
            if not df.empty:
                model_predicted_csv = df.to_csv(index=False)
                process_csv(model_predicted_csv)
        except Exception as e:
            print('model error!')
    else:
        df = pd.DataFrame()  # Create an empty DataFrame if no data is found

    # Convert combined data to a DataFrame
    combined_df = pd.DataFrame(combined_data)
    combined_df.drop_duplicates(subset=['refseq_ids', 'variation_ids'], inplace=True)

    # Combine the combined_df and df
    if not df.empty:
         combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Create a CSV file for combined data
    combined_csv_buffer = io.StringIO()
    combined_df.to_csv(combined_csv_buffer, index=False)
    combined_csv_buffer.seek(0)

    # Send the email with the CSV attachment
    try:
        email_message = EmailMessage(
            'Your files are here',
            'Hello, here are your files.',
            to=[email]
        )
        email_message.attach('combined_data.csv', combined_csv_buffer.getvalue(), 'text/csv')
        email_message.send()
        return 'success'
    except Exception as e:
        return f'Error: {str(e)}'


def upload(request):
    if request.method == 'POST':
        files = request.FILES.getlist('files')
        for f in files:
            process_csv(f)
        return render(request, 'upload.html', {'msg': "Upload Successfully !"})

    return render(request, 'upload.html')


def process_csv(csv_input, batch_size=10000):
    "process_csv Function"
    print("process function")

    # Determine if csv_input is a string or a file-like object
    if isinstance(csv_input, str):
        # If it's a string, wrap it in a StringIO object
        csv_input = io.StringIO(csv_input)
    elif hasattr(csv_input, 'read'):
        # If it's a file-like object, reset to the start if needed
        csv_input.seek(0)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_input)
    print("in process function", df.head())

    # List to hold all the objects to be created
    objects_to_create = []

    # Iterate over the rows of the DataFrame and prepare to save to the database if not duplicate
    for _, row in df.iterrows():
        # Check for duplication based on refseq_ids and variation_ids
        if not data.objects.filter(refseq_ids=row[0], variation_ids=row[1]).exists():
            objects_to_create.append(data(
                refseq_ids=row[0],
                variation_ids=row[1],
                meanProb=row[2],
                stdProb=row[3],
                pred_label=row[4]
            ))

        # Bulk create objects in batches
        if len(objects_to_create) >= batch_size:
            with transaction.atomic():
                data.objects.bulk_create(objects_to_create)
            objects_to_create = []

    # Create any remaining objects
    if objects_to_create:
        with transaction.atomic():
            data.objects.bulk_create(objects_to_create)

def about(request):
    return render(request, 'about.html')


def disclaimer(request):
    return render(request, 'disclaimer.html')


def loading_dictionary():
    genes_dict_path = os.path.join(settings.DATA_FILES_DIR, 'Updated_Ens_genes_dict.txt')
    protein_gene_info = {}
    with open(genes_dict_path, 'r') as file:
        for line in file:
            line = line.strip()
            parts = line.split('\t')
            if len(parts) == 2:
                prot_id, gene_name = parts[0], parts[1]
                protein_gene_info[prot_id] = gene_name
    print("Dictionary loaded successfully.")

    go_file_path = os.path.join(settings.DATA_FILES_DIR, 'gene_GO_sequence.csv')
    go_file = pd.read_csv(go_file_path)
    go_ratio_dict = dict(zip(go_file['ensembl_gene_ids'], go_file['GO_ratio']))

    select_indices_file_path = 'selected_feat_indecs_50.csv'

    models_dir = settings.MODELS_DIR

    return go_ratio_dict, models_dir, select_indices_file_path, protein_gene_info


def get_prediction(record_file):
     # = pd.read_csv(file_name)
    go_ratio_dict, models_dir, select_indices_file_path, protein_gene_info = loading_dictionary()

    test_feature = extract_feature(record_file)
    select_indices_file = pd.read_csv("web/selected_feat_indecs_50.csv")
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
    print('after for loop line')
    go_data = pd.DataFrame({'refseq_ids': refseq_ids, 'GO_ratio': go_ratio_value})
    test_data['GO_value'] = pd.to_numeric(go_data['GO_ratio'])
    test_data_wGO = test_data.values

    n_iterations = 200
    pred_matrix = np.empty((len(test_data_wGO), 0))

    for i in range(n_iterations):
        model_filename = os.path.join("web/models/", f'lgbm_model_{i + 1}.joblib')
        lgbm = joblib.load(model_filename)
        pred_probs = lgbm.predict_proba(test_data_wGO)[:, 1]
        pred_matrix = np.column_stack([pred_matrix, pred_probs])

    mean_probs = np.mean(pred_matrix, axis=1)
    pred_labels = np.where(mean_probs >= 0.5, 'P', 'N')
    stdprob = np.std(pred_matrix, axis=1)
    print('after for loop 2 line')
    results_table = pd.DataFrame({
        "refseq_ids": record_file["refseq_ids"],
        "variation_ids": record_file["variation_ids"],
        "meanProb": np.round(mean_probs, 3),
        "StandardDev": np.round(stdprob, 3),
        "pred_class": pred_labels,
    })

    table_headers = ["refseq_ids", "variation_ids", "meanProb", "StandardDev", "pred_class"]
    table_data = results_table.values.tolist()

    # Print results in table format for debugging
    # print(tabulate(table_data, headers=table_headers, floatfmt=".3f"))
    print('printing before returning')
    return results_table


def Get_pred_with_model(email, file):
    combined_df = get_prediction(file)
    csv_buffer = io.StringIO()
    combined_df.to_csv(csv_buffer, index=False)
    csv_buffer.seek(0)  # Go to the start of the StringIO buffer

    # Send the email with CSV attachment
    try:
        email_message = EmailMessage(
            'Your files are here',
            'Hello, this is your file',
            to=[email]
        )
        email_message.attach('combined_data.csv', csv_buffer.getvalue(), 'text/csv')
        email_message.send()
        return 'success'
    except Exception as e:
        print(f"Error: {e}")
        return 'Error'
