import sqlite3
import csv
import glob
import os
import random
import glob
import src.data_handling as dh
import pandas as pd
import secrets
import string
import markdown
import seaborn as sns
import matplotlib.pyplot as plt


def find_sqlite(input_dir):
    return glob.glob(os.path.join(input_dir, "*.sqlite"))


def xtract_from_sql(db_name, path_to_sqlite):
    # Check if db_name starts with a number, if so, add 'sql_' to the start
    db_name = f"sql_{db_name}"
    db_name = db_name.replace('-', '_')

    # Connect to the SQLite database
    conn = sqlite3.connect(path_to_sqlite)
    cursor = conn.cursor()
    cursor.execute(f"DROP TABLE IF EXISTS {db_name}_extracted;")

    # Create a new table with the required columns
    cursor.execute(f"""
    CREATE TABLE IF NOT EXISTS {db_name}_extracted (
        attributes_accession TEXT,
        neighbors_accession TEXT,
        attributes_start INTEGER,
        attributes_stop INTEGER,
        seq_len INTEGER,
        cluster_num INTEGER,
        organism TEXT,
        family_desc TEXT,
        neighbors_start INTEGER,
        neighbors_stop INTEGER,
        a_direction TEXT,
        n_direction TEXT,
        nseq_len INTEGER,
        rel_start INTEGER,
        rel_stop INTEGER,
        a_family TEXT
    )
    """)

    # Insert required data from the 'attributes' and 'neighbors' tables
    cursor.execute(f"""
    INSERT INTO {db_name}_extracted (attributes_accession, neighbors_accession, attributes_start, attributes_stop, seq_len, cluster_num, organism, family_desc, neighbors_start, neighbors_stop, a_direction, n_direction, nseq_len, rel_start, rel_stop, a_family)
    SELECT a.accession, n.accession, a.start, a.stop, a.seq_len, a.cluster_num, a.organism, n.ipro_family_desc, n.start, n.stop, a.direction, n.direction, n.seq_len, n.rel_start, n.rel_stop, a.ipro_family
    FROM attributes a
    JOIN neighbors n
    ON a.sort_key = n.gene_key
    """)

    # Commit the changes
    conn.commit()

    # Fetch the data from the new table
    data_df = pd.read_sql_query(
        "SELECT * FROM " + db_name + "_extracted", conn)

    # Close the connection
    conn.close()

    return data_df


def extract_sqlite_paths(input_folder, output_folder):
    sql_files = find_sqlite(input_folder)
    csv_names = []
    output_folders = []
    fs = []
    db_names = []
    for f in sql_files:
        db_name = os.path.basename(f).split(".", 1)[0]
        output_folder = dh.init_path(output_folder, db_name)
        csv_name = os.path.join(output_folder, f"{db_name}_extracted")
        # csv_name = xtract_from_sql(output_folder, db_name, f)
        fs.append(f)
        db_names.append(db_name)
        csv_names.append(csv_name)
        output_folders.append(output_folder)
    return output_folders, csv_names, fs, db_names


def convert_sqlite_data(df, output_folder, csv_name, node_csv_folder, n_n=4, from_stop=True):
    '''CSV should be sql_*_extracted.csv. User can specify how many upstream and downstream genes to include from flask app.'''
    node_csv_path = glob.glob(os.path.join(
        node_csv_folder, "*node_table.csv"))
    print(node_csv_path)
    if node_csv_path:
        ndf = pd.read_csv(node_csv_path[0])
        ndf = ndf[["List of IDs in Rep Node",
                   "Phylum", "Order", "InterPro (Domain)"]]
        # Split the 'list_of_ids' column on '|', stack the resulting DataFrame, and reset the index
        ndf_ids = ndf["List of IDs in Rep Node"].str.split("|", expand=True).stack(
        ).reset_index(level=1, drop=True).to_frame('list_of_ids')
        # Reset the index of the original DataFrame
        ndf.reset_index(inplace=True, drop=True)
        # Join the expanded 'list_of_ids' column with the other columns in the original DataFrame
        node_df = ndf_ids.join(ndf.drop('List of IDs in Rep Node', axis=1))
        df = df.merge(node_df, left_on='attributes_accession',
                      right_on='list_of_ids')
        df = df[~df['Order'].str.contains('\|')]

    df['same_direction'] = df['a_direction'] == df['n_direction']
    # by default (from_stop = True) measure closest neighbor by relative stop position. Better for capturing upstream genes
    df['abs_rel_start'] = df['rel_start'].abs()
    df['abs_rel_stop'] = df['rel_stop'].abs()
    if from_stop:
        df_sorted = df.sort_values(['attributes_accession', 'abs_rel_stop'])
    else:
        df_sorted = df.sort_values(['attributes_accession', 'abs_rel_start'])

    df_op = df_sorted[df_sorted['same_direction'] == True]

    # Group by attributes_accession and take the first 4 rows of each group
    df_filtered = df_op.groupby('attributes_accession').head(n_n-1)

    # Group the DataFrame by attributes_accession and aggregate neighbors_accession and family_desc
    grouped = df_filtered.groupby('attributes_accession').agg({
        'neighbors_accession': list,
        'family_desc': list,
        'nseq_len': list
    }).reset_index()

    temp_df = pd.DataFrame(
        grouped['neighbors_accession'].tolist(), index=grouped.index)
    num_columns = temp_df.shape[1]
    grouped[[f'neighbors_accession_{i}' for i in range(
        1, num_columns + 1)]] = temp_df
    temp_df_family_desc = pd.DataFrame(
        grouped['family_desc'].tolist(), index=grouped.index)
    num_columns_family_desc = temp_df_family_desc.shape[1]
    grouped[[f'family_desc_{i}' for i in range(
        1, num_columns_family_desc + 1)]] = temp_df_family_desc
    temp_df_nseq_len = pd.DataFrame(
        grouped['nseq_len'].tolist(), index=grouped.index)
    num_columns_nseq_len = temp_df_nseq_len.shape[1]
    grouped[[f'nseq_len_{i}' for i in range(
        1, num_columns_nseq_len + 1)]] = temp_df_nseq_len

    # Drop the original columns
    grouped = grouped.drop(
        columns=['neighbors_accession', 'family_desc', 'nseq_len'])

    # Merge the organism column and cluster_num back into the final table
    final_table = grouped.merge(
        df[['attributes_accession', 'organism',
            'cluster_num', 'seq_len', 'Phylum', 'Order']].drop_duplicates(),
        on='attributes_accession'
    )

    # Function to generate a random 6-character string (letters and numbers)
    def generate_group_code():
        characters = string.ascii_letters + string.digits
        return ''.join(secrets.choice(characters) for _ in range(6))

    # Add a new column 'group_code' with a randomly generated 6-character string
    final_table['group_code'] = final_table.apply(
        lambda _: generate_group_code(), axis=1)
    final_table['Genus'] = final_table['organism'].str.split(" ").str[0]
    # csv_name = os.path.basename(csv_path).split(".", 1)[0]
    # csv_out = os.path.join(f"{csv_name}_closest_neighbors.csv")
    # # Save the result to a new CSV file
    # final_table.to_csv(csv_out, index=False)

    return final_table


def convert_sqlite(input_folder, output_folder, n_n=4):
    # output folder is output/extracted_sqlite
    sql_files = dh.list_csvs(output_folder)
    new_csv_paths = []
    for f in sql_files:
        # print(f)
        # print("!!!\n!!!!!!!!!!!!\nxxxsddfsfdsf\n!!!!!!!!!!!!\n!!!")
        db_name = os.path.basename(f).split(".", 1)[0]
        new_csv_path = convert_sqlite_data(
            f, input_folder, db_name, n_n)
        new_csv_paths.append(new_csv_path)
    return new_csv_paths


def sample_and_search(df, csv_name, output_folder, input_folder, accessions_folder, num_n=1, min_cluster_size=30, prop_size=1000):
    '''
    Big mess. Point is to convert the big neighbors csv file into a reduced dataset

    This randomly samples the nieghbours.csv file proportionally to the number of rows belonging to each gene cluster. However, the minimum number of smaples is 10, so the data will be a bit biased by these.  However the far away clsuters can simply be thrown out if they ruin an alignment. This is so that alignments/some inference can be made with all the clusters present in the data. 

    '''
    # Tabulate the count of each unique cluster number
    cluster_counts = df['cluster_num'].value_counts()
    count_total = cluster_counts.sum()
    # print(f"Clusters: \n {cluster_counts.head(15)}")
    # normalize dataset to 1000 rows by default
    norm_factor = prop_size / count_total

    # Sample rows by cluster, proportional to the number of rows in each cluster
    sampled_rows = []
    for cluster_num, count in cluster_counts.items():
        # If the number of rows in a cluster is less than min (set in function header), don't sample

        if count < min_cluster_size:
            continue

        # Calculate the sample size based on the normalized count
        # You can adjust the factor (e.g., 100 for proportions to represent %) to control the sample size
        sample_size = int((norm_factor) * count)
        # Prevent tiny useless datasets
        if sample_size < 10:
            sample_size = 10
        # Ensure that the sample size is at least 1
        sample_size = max(1, sample_size)
        print(
            f'count: {count}, count_total: {count_total}, sample_size: {sample_size}')
        # Sample rows from the current cluster
        cluster_sample = df[df['cluster_num']
                            == cluster_num].sample(sample_size)

        # Append the sampled rows to the list
        sampled_rows.extend(cluster_sample.to_dict("records"))
    sampled_df = pd.DataFrame(sampled_rows)

    # Also randomly sample 10 rows
    random_sample_rows = []
    small_random_sample_rows = []
    for cluster_num, count in cluster_counts.items():
        # Check if the cluster exists in the dataframe
        if count < min_cluster_size:
            continue
        # Sample 10 and 2 random rows from the current cluster
        cluster_random_sample = df[df['cluster_num'] == cluster_num].sample(10)
        cluster_random_sample2 = cluster_random_sample.sample(1)
        # Append the sampled rows to the list
        random_sample_rows.extend(cluster_random_sample.to_dict("records"))
        small_random_sample_rows.extend(
            cluster_random_sample2.to_dict("records"))

    # Save the samples to output named {csv_name}_cluster{cluster_num}_t10rsample.csv
    random_sample_df = pd.DataFrame(random_sample_rows)
    small_random_sample_df = pd.DataFrame(small_random_sample_rows)
    # random_sample_df.to_csv(os.path.join(
    #     output_folder, f"{csv_name}_cluster_t10rsample.csv"), index=False)

    # Search for accessions txt file
    accessions_files = glob.glob(os.path.join(
        accessions_folder, "*accessions.txt"))

    combined_df = pd.DataFrame()
    combined_df2 = pd.DataFrame()
    combined_df3 = pd.DataFrame()
    for accessions_file in accessions_files:
        print("Found accessions.txt file")
        with open(accessions_file, "r") as f:
            accessions = [line.strip() for line in f.readlines()]

        # Match the accessions to the 'attributes_accession' column and extract the rows
        selected_rows = df[df['attributes_accession'].isin(accessions)]
        combined_df = pd.concat(
            [sampled_df, selected_rows], ignore_index=True).drop_duplicates()

        combined_df2 = pd.concat(
            [random_sample_df, selected_rows], ignore_index=True).drop_duplicates()

        combined_df3 = pd.concat(
            [small_random_sample_df, selected_rows], ignore_index=True).drop_duplicates()
        # Save the matching rows as {csv_name}_selected.csv
        # csv_name = os.path.basename(input_csv_path).split(".", 1)[0]
        selected_rows.to_csv(f"{csv_name}_selected.csv", index=False)

        # Find accessions not found in the input csv
        not_found_accessions = set(
            accessions) - set(selected_rows['attributes_accession'].unique())

        # Save the not found accessions to {csv_name}_NotFound.txt
        if not_found_accessions:
            with open(f"{csv_name}_NotFound.txt", "w") as nf:
                nf.write("\n".join(not_found_accessions))

    # If no accession text file in input folder, will just take the samples
    if combined_df.empty:
        combined_df = sampled_df
    if combined_df2.empty:
        combined_df2 = random_sample_df
    if combined_df3.empty:
        combined_df3 = small_random_sample_df
    # Write out dataframes
    combined_df.to_csv(f"{csv_name}_proportional_rsample.csv", index=False)
    combined_df2.to_csv(f"{csv_name}_big_rsample.csv", index=False)
    combined_df3.to_csv(f"{csv_name}_small_rsample.csv", index=False)
    # Output smaller list of accessions for Program 2 to run on
    filtered_df3 = combined_df3.filter(regex='_accession')
    filtered_df3.to_csv(
        f"{csv_name}_small_sample_neighbor_accessions.csv", index=False)

    with open(f"{csv_name}_proportionial_query_accs.txt", "w") as nf:
        nf.write("\n".join(combined_df['attributes_accession'].unique()))
    # return the proportional dataframe to be matched with uniprot data.
    return combined_df


def main():
    print("No test implemented right now")
    return


if __name__ == "__main__":
    main()
