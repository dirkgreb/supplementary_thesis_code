import subprocess
from Bio import AlignIO
from Bio.Align import AlignInfo

def calculate_conservation(alignment_file, output_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)

    with open(output_file, "w") as output:
        for i in range(len(alignment[0])):
            position = summary_align.pos_specific_score(i)
            first_seq_aa = alignment[0][i]
            score = position.get(first_seq_aa, 0)
            output.write(f"{i + 1}\t{first_seq_aa}\t{score}\n")

def run_muscle(input_file, output_file, maxiters=16, diags=True):
    muscle_cmd = [
        "muscle",
        "-in", input_file,
        "-out", output_file,
        "-maxiters", str(maxiters),
    ]

    if diags:
        muscle_cmd.append("-diags")

    subprocess.run(muscle_cmd, check=True)


# Example usage
input_fasta = "input.fasta"  # Path to your input fasta file
aligned_sequences = "aligned.aln"  # Path to save the aligned sequences
conservation_scores = "conservation_scores.txt"  # Path to save conservation scores

run_muscle(input_fasta, aligned_sequences)
calculate_conservation(aligned_sequences, conservation_scores)
