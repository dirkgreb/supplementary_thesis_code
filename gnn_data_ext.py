import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import src.data_handling as dh
import src.run_get_accessions as rga
import src.concurrent_system_calls as csc
import glob
import json
import requests
from Bio import Entrez, SeqIO

# Always provide your email when using NCBI Entrez API
Entrez.email = "dirkgrebenc@gmail.com"


def parse_uniprot_record(record):
    if len(record) == 0:
        return None
    record = record[0]  # Access the first element in the list
    accession = record.get("accession", "Not found")
    gene = record.get("gene", [{}])[0].get(
        "orfNames", [{}])[0].get("value", "Not found")

    organism = record.get("organism", {})
    organism_name = organism.get("names", [{}])[0].get("value", "Not found")
    lineage = organism.get("lineage", [])

    protein = record.get("protein", {}).get("submittedName", [{}])[
        0].get("fullName", {}).get("value", "Not found")

    db_references = record.get("dbReferences", [])
    caution = record.get("caution", [{}])[0].get("value", "Not found")
    gene_id = "Not found"
    RefSeq = "Not found"
    AlphaFoldDB = "Not found"
    embl = "Not Found"
    for ref in db_references:
        # print(ref["type"])
        if ref["type"] == "GeneID":
            gene_id = ref["id"]
        if ref["type"] == "RefSeq":
            RefSeq = ref["id"]
        if ref["type"] == "AlphaFoldDB":
            AlphaFoldDB = ref["id"]
        if ref["type"] == "EMBL":
            embl = ref["id"]

    sequence = record.get("sequence", {}).get("sequence", "Not found")

    return {
        "accession": accession,
        "gene_name": gene,
        "organism_name": organism_name,
        "lineage": lineage,
        "protein_name": protein,
        "gene_id": gene_id,
        "refseq": RefSeq,
        "alphafold": AlphaFoldDB,
        "embl": embl,
        "caution": caution,
        "sequence": sequence
    }


def get_uniprot_record(accession_codes: list):
    print(f"Fetching data for accession codes: {accession_codes} ")
    data = []

    for code in accession_codes:
        url = f"https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession={code}"

        response = requests.get(url, headers={"Accept": "application/json"})

        if response.status_code == 200:
            record = response.json()
            parsed_data = parse_uniprot_record(record)  # function from above
            if parsed_data is not None:
                data.append(parsed_data)

        else:
            print(
                f"Failed to fetch data for accession code {code}, status code: {response.status_code}")

    if data:
        df = pd.DataFrame(data)
        print("Done!")
    else:
        df = pd.DataFrame()
        print("Data didn't go into DataFrame")

    return df


def extract_accessions_from_prop(input_paths: list):
    return dh.convert_n_csvs_to_lists(input_paths)


def accs_dfs_to_uniprot(dfs, csv_name="*_selected.csv"):
    for df in dfs:
        # cleans up extra info from prop csv

        # print(trunc_rows)
        total = len(trunc_rows)
        done = 1
        datas = []
        for row in trunc_rows:
            print(f"Getting records for {row}")
            data = get_uniprot_record(row)
            if not data.empty:
                datas.append(data)
            print(
                f"Got it! \n {done} of {total} records -> {round(done/total*100, 2)}%\n")
            done += 1
        result = pd.concat(datas, ignore_index=True)
        dfs.append(result)
        csv_n = csv_name[1:]
        output = os.path.join(input_f, f"{f}_uniprot_data.csv")
        result.to_csv(os.path.join(f"{output}"), index=False)
    return dfs


def esql_to_uniprot(input_folder='output/sqlite_data', csv_name="*_selected.csv"):
    folders = [x for x in dh.directory_list(input_folder)]
    dfs = []
    for f in folders:
        dh.list_csvs(os.path.join(input_folder, f))
        input_f = os.path.join(input_folder, f)
        props = glob.glob(os.path.join(input_f, csv_name))
        rows = extract_accessions_from_prop(props)
        # cleans up extra info from prop csv
        trunc_rows = [x[0:4] for x in rows[0]]
        # print(trunc_rows)
        total = len(trunc_rows)
        done = 1
        datas = []
        for row in trunc_rows:
            print(f"Getting records for {row}")
            data = get_uniprot_record(row)
            if not data.empty:
                datas.append(data)
            print(
                f"Got it! \n {done} of {total} records -> {round(done/total*100, 2)}%\n")
            done += 1
        result = pd.concat(datas, ignore_index=True)
        dfs.append(result)
        csv_n = csv_name[1:]
        output = os.path.join(input_f, f"{f}_uniprot_data.csv")
        result.to_csv(os.path.join(f"{output}"), index=False)
    return dfs


def sdf_to_uniprot(df, n_n):
    fdf = df.filter(regex="_accession", axis=1)
    code_rows = dh.rows_to_lists(fdf)
    datas = [get_uniprot_record(code_row) for code_row in code_rows]
    udf = pd.concat(datas, ignore_index=True)
    table = combine_sample_uniprot_dfs(df, udf, n_n)
    return table


def combine_sample_uniprot_dfs(search_sample, udf, max_neighbors):

    ex_df = udf[['accession', 'lineage', 'sequence', 'gene_id', 'caution']]
    ex_df = ex_df.copy()
    ex_df[['kingdom', 'phylum', 'class', 'order', 'family', 'genus']] = pd.DataFrame(
        ex_df['lineage'].astype(str).apply(eval).tolist(), index=ex_df.index).iloc[:, :6]
    ex_df.drop(columns=['lineage'], inplace=True)

    # Extract from out the neighbor accessions, rename columns and then recombine
    df_list = []
    df0 = search_sample[['attributes_accession', 'seq_len', 'organism', 'cluster_num', 'group_code']].rename(
        columns={'attributes_accession': 'accession'})
    df0['family'] = 'query'
    df0['position'] = 0
    df_list.append(df0)

    # Update this value with the maximum number of neighbors you want to process
    for i in range(1, max_neighbors):
        df_temp = search_sample[
            [f'neighbors_accession_{i}', f'family_desc_{i}', f'nseq_len_{i}', 'organism', 'cluster_num', 'group_code']].rename(
            columns={f'neighbors_accession_{i}': 'accession', f'family_desc_{i}': 'family', f'nseq_len_{i}': 'seq_len'})
        df_temp['position'] = i
        df_list.append(df_temp)

    gn_df = pd.concat(df_list, axis=0)

    # Join the ex_df and gn_df
    final_table = gn_df.merge(ex_df, on='accession')
    return final_table


def combine_uniprot_sample_tables(sql_folder_name, sql_folder):
    # DEPRECIATED
    '''
    Combine sqlite dataframes into one dataframe
    sql_folder_name: name of the folder that contains the csvs of interest
    sql_folder: full path to the folder that contains the csvs of interest
    '''

    # Find the files, convert to dataframes
    search_sample_path = glob.glob(os.path.join(
        sql_folder, "*_search_sample.csv"))[0]
    uniprot_data_path = glob.glob(os.path.join(
        sql_folder, "*_uniprot_data.csv"))[0]

    search_sample = pd.read_csv(search_sample_path)
    uniprot_data = pd.read_csv(uniprot_data_path)

    # Extract from uniprot_data
    ex_df = uniprot_data[['accession', 'lineage', 'sequence']]
    ex_df = ex_df.copy()
    ex_df[['kingdom', 'phylum', 'class', 'order', 'family', 'genus']] = pd.DataFrame(
        ex_df['lineage'].apply(eval).tolist(), index=ex_df.index).iloc[:, :6]
    ex_df.drop(columns=['lineage'], inplace=True)

    # Extract from out the neighbor accessions, renemae columns and then recombine
    df0 = search_sample[['attributes_accession', 'seq_len', 'organism',
                         'cluster_num', 'group_code']].rename(
                             columns={'attributes_accession': 'accession'})
    df1 = search_sample[['neighbors_accession_1', 'family_desc_1', 'nseq_len_1', 'organism', 'cluster_num', 'group_code']].rename(
        columns={'neighbors_accession_1': 'accession', 'family_desc_1': 'family', 'nseq_len_1': 'seq_len'})
    df2 = search_sample[['neighbors_accession_2', 'family_desc_2', 'nseq_len_2', 'organism', 'cluster_num', 'group_code']].rename(
        columns={'neighbors_accession_2': 'accession', 'family_desc_2': 'family', 'nseq_len_2': 'seq_len'})
    df3 = search_sample[['neighbors_accession_3', 'family_desc_3', 'nseq_len_3', 'organism', 'cluster_num', 'group_code']].rename(
        columns={'neighbors_accession_3': 'accession', 'family_desc_3': 'family', 'nseq_len_3': 'seq_len'})
    df0['family'] = 'query'
    df0['position'] = 0
    df1['position'] = 1
    df2['position'] = 2
    df3['position'] = 3

    gn_df = pd.concat([df0, df1, df2, df3], axis=0)

    # Join the ex_df and gn_df
    final_table = gn_df.merge(ex_df, on='accession')
    final_table.to_csv(os.path.join(
        sql_folder, f"{sql_folder_name}_ds_master.csv"), index=False)

    return final_table


def main():
    print("No tests yet")
    return


if __name__ == "__main__":
    main()
