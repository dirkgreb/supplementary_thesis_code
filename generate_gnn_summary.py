import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import markdown


def summarize_clusters(df):
    cluster_counts_full = df['cluster_num'].value_counts().rename_axis(
        'Cluster_Number').reset_index(name='counts')
    norm_counts = df['cluster_num'].value_counts(normalize=True).rename_axis(
        'Cluster_Number').reset_index(name='norm_counts').head(15)
    cluster_counts = cluster_counts_full.head(15)
    both = cluster_counts.merge(norm_counts, on='Cluster_Number')
    count_total = cluster_counts_full.sum()
    summary_string = f"Cluster Counts\n"
    summary_string += f"\n| Cluster | Count | Percentage |\n"
    summary_string += f"| ---- | ---- | ---- |\n"
    for _, row in both.iterrows():
        summary_string += f"| {row['Cluster_Number']} | {row['counts']} | {round(row['norm_counts']*100, 2)}% |\n"
    # summary_string += f"15 largest clusters: \n {cluster_counts.head(15)}\n\n"
    summary_string += f"\n"
    return summary_string


def summarize_neighbor_fams(df, col_prefix='family_desc'):
    # Summarize the most frequent words in the family_desc columns
    family_desc_columns = [
        col for col in df.columns if col.startswith(col_prefix)]
    summary_string = ""
    # Process each family_desc column
    for col in family_desc_columns:
        # Split values in the column by ';'
        # pos = df[col].str.split(';', expand=True).dropna()
        # pos_vals = pos.stack().value_counts()

        pos = df[col]
        pos_vals = pos.value_counts()

        # Calculate the total number of unique values
        unique_values_count = len(pos_vals)

        # Get the top 10 most common values
        top_ten_values = pos_vals.head(10)
        top_ten_values = df[col].value_counts().rename_axis(
            f'{col}').reset_index(name='counts').head(10)
        top_ten_values_n = df[col].value_counts(normalize=True).rename_axis(
            f'{col}').reset_index(name='counts_n').head(10)
        top_ten_values = top_ten_values.merge(top_ten_values_n, on=col)

        # save the summary for the current column
        summary_string += str(
            f"{col.upper()} \n Total number of unique values: {unique_values_count}  \n")
        summary_string += str(
            f"Top ten most common {col} and their frequencies:\n")
        summary_string += str(f"\n| {col} | Count | Percent |\n")
        summary_string += str(f"| ---- | ---- | ---- |\n")

        for _, row in top_ten_values.iterrows():
            if row[col] == ' ' or row[col] == '' or row[col] == 'nan' or row[col] == 'None':
                row[col] = 'No PFAM Annotation'
            summary_string += str(
                f"| {row[col]} | {int(row['counts'])} | {round(row['counts_n']*100, 2)}% |\n")
        # for word, count in top_ten_values.items():
        #     out_string += str(f"| {word} | {count} |\n")
        summary_string += str("\n\n")

        # summary_string += str(f"\n| {col} | Frequency |\n")
        # summary_string += str(f"| ---- | ---- |\n")
        # for word, count in top_ten_values.items():
        #     if word == ' ' or word == '' or word == 'nan' or word == 'None':
        #         word = 'No PFAM Annotation'
        #     summary_string += f"| {word} | {count} |\n"
        # summary_string += f"\n\n"

    return summary_string


def summarize_phylo(df, col_prefix):
    # Summarize the most frequent words in the family_desc columns
    phylo_columns = [col for col in df.columns if col == col_prefix]
    # Process each family_desc column
    out_string = ""
    for col in phylo_columns:
        # Calculate the total number of unique values
        unique_values_count = len(df[col].value_counts())

        # Get the top ten most common values
        top_ten_values = df[col].value_counts().rename_axis(
            f'{col}').reset_index(name='counts').head(10)
        top_ten_values_n = df[col].value_counts(normalize=True).rename_axis(
            f'{col}').reset_index(name='counts_n').head(10)
        top_ten_values = top_ten_values.merge(top_ten_values_n, on=col)
        # Print the summary for the current column
        out_string += str(
            f"{col.upper()} \n Total number of unique values: {unique_values_count}  \n")
        out_string += str(
            f"Top ten most common {col} and their frequencies:\n")
        out_string += str(f"\n| {col} | Count | Percent |\n")
        out_string += str(f"| ---- | ---- | ---- |\n")

        for _, row in top_ten_values.iterrows():
            out_string += str(
                f"| {row[col]} | {int(row['counts'])} | {round(row['counts_n']*100, 2)}% |\n")
        # for word, count in top_ten_values.items():
        #     out_string += str(f"| {word} | {count} |\n")
        out_string += str("\n\n")
    return out_string


def save_pie_chart(df, column, output_path, cluster_num, threshold=0.02):
    # extract the data for the specified column
    column_data = df[column]
    value_counts = column_data.value_counts(normalize=True)

    # Group values below the threshold into an "Other" category
    mask = value_counts >= threshold
    other_value = value_counts[~mask].sum()
    value_counts = value_counts[mask]
    if other_value > 0:
        value_counts["Other"] = other_value

    values = value_counts.tolist()
    labels = value_counts.index.tolist()

    # generate a color palette using Seaborn
    palette = sns.color_palette("Set2", len(labels))

    # convert the palette to a dictionary that maps labels to colors
    color_dict = dict(zip(labels, palette))

    # extract the colors for each label
    colors = [color_dict[label] for label in labels]

    # create the pie chart with percentage labels
    plt.pie(values, labels=labels, colors=colors, autopct='%1.1f%%')

    # add a title
    plt.title(f"{column} Cluster {cluster_num} Distribution")
    output_name = os.path.join(
        output_path, f'cluster_{cluster_num}_{column}_pie_chart.png')
    # save the plot to a file
    plt.savefig(output_name)
    plt.close()
    figure_name = os.path.basename(output_name)
    figure_folder = os.path.basename(os.path.dirname(output_name))
    figure_path = os.path.join(figure_folder, figure_name)
    # figure_path = output_name.split('/', 2)[2]
    # return str(f'./{figure_path}')
    return str(f'{figure_path}')


def save_numerical_hist(df, column, output_path, cluster_num):
    plt.figure()
    sns.histplot(df[column].dropna(), kde=False)
    plt.title(f'Histogram of {column} for Cluster {cluster_num}')
    plt.xlabel(column)
    plt.ylabel('Frequency')
    output_name = os.path.join(
        output_path, f'cluster_{cluster_num}_{column}_hist.png')
    plt.savefig(output_name)
    plt.close()
    figure_name = os.path.basename(output_name)
    figure_folder = os.path.basename(os.path.dirname(output_name))
    figure_path = os.path.join(figure_folder, figure_name)
    # return str(f'./{figure_path}')
    return str(f'{figure_path}')


def generate_report(df, csv_name, report_out_path):
    os.makedirs(report_out_path, exist_ok=True)
    figures_path = os.path.join(report_out_path, "figures")
    os.makedirs(figures_path, exist_ok=True)
    cluster_nums = sorted(df['cluster_num'].unique())
    numerical_columns = [col for col in df.columns if "seq_len" in col]

    out_text = []
    out_text.append(f'# {" ".join(csv_name.split("_"))} Report\n\n')
    out_text.append(
        f'## Overall Dataset Stats \n\n ------------------- \n\n Total number of operons: {len(df)}\n\n')
    out_text.append(summarize_clusters(df))
    out_text.append(summarize_phylo(df, 'Phylum'))
    phyl_chart = save_pie_chart(df, 'Phylum', figures_path, 'all')
    out_text.append(f'![Phylum Pie Chart for All Clusters]({phyl_chart})\n\n')
    out_text.append(summarize_phylo(df, 'Order'))
    order_chart = save_pie_chart(df, 'Order', figures_path, 'all')
    out_text.append(f'![Order Pie Chart for All Clusters]({order_chart})\n\n')

    out_text.append(summarize_phylo(df, 'Genus'))
    genus_chart = save_pie_chart(df, 'Genus', figures_path, 'all')
    out_text.append(f'![Order Pie Chart for All Clusters]({genus_chart})\n\n')
    # out_text.append(summarize_phylo(df, 'Genus'))
    out_text.append(summarize_neighbor_fams(df, 'family_desc_'))
    out_text.append(f'## By Cluster Stats for clusters with >20 members  \n\n')
    for cluster_num in cluster_nums:
        cdf = df[df['cluster_num'] == cluster_num]
        if len(cdf) < 20:
            continue

        # Save pie charts and histogram for the current cluster
        phylo_chart_p = save_pie_chart(
            cdf, 'Phylum', figures_path, cluster_num)
        order_chart_p = save_pie_chart(cdf, 'Order', figures_path, cluster_num)
        genus_chart_p = save_pie_chart(cdf, 'Genus', figures_path, cluster_num)
        histogram_paths = [save_numerical_hist(
            cdf, num_col, figures_path, cluster_num) for num_col in numerical_columns]

        out_text.append(f'### Cluster {cluster_num}  \n\n')
        out_text.append(f'#### Phylogenic Summary  \n\n')
        out_text.append(str(summarize_phylo(cdf, 'Phylum')))
        out_text.append(
            f'![Phylum Pie Chart for Cluster {cluster_num}]({phylo_chart_p})\n\n')
        out_text.append(str(summarize_phylo(cdf, 'Order')))
        out_text.append(
            f'![Order Pie Chart for Cluster {cluster_num}]({order_chart_p})\n\n')
        out_text.append(str(summarize_phylo(cdf, 'Genus')))
        out_text.append(
            f'![Genus Pie Chart for Cluster {cluster_num}]({genus_chart_p})\n\n')
        out_text.append(f'#### Neighbouring Genes Summary\n\n')
        out_text.append(str(summarize_neighbor_fams(
            cdf, col_prefix="family_desc_")))
        # cdf['unique_families'] = cdf.apply(lambda fam_desc: cdf[f'{fam_desc}'].unique(
        # ) for fam_desc in cdf.columns if fam_desc.startswith("family_desc"))
        
        for his_path in histogram_paths:
            out_text.append(
                f'![Seq_len_hist Histogram for Cluster {cluster_num}]({his_path})\n\n')
    return out_text


def write_report(db_names, sql_outs, tables, csv_names, n_output_folder):
    report_paths = [os.path.join(
        n_output_folder, f'{db_name}_report.md') for db_name, sql_out in zip(db_names, sql_outs)]
    reports = (generate_report(table, csv_name, n_output_folder)
               for table, csv_name in zip(tables, csv_names))
    for report, rep_path in zip(reports, report_paths):
        with open(rep_path, 'w') as f:
            f.writelines(map(str, report))
    return report_paths


def main():
    print("not tests yet")
    return


if __name__ == "__main__":
    main()
