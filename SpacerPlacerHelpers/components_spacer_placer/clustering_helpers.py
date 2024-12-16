import csv
from collections import defaultdict
import Levenshtein


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence, handling ambiguous bases."""
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N'
    }
    return ''.join(complement[base] for base in reversed(sequence))


def is_edit_distance_one_or_zero(seq1, seq2):
    """Checks if the edit distance between two sequences is 0 or 1."""
    if abs(len(seq1) - len(seq2)) > 1:
        return False
    direct_distance = Levenshtein.distance(seq1, seq2)
    reverse_distance = Levenshtein.distance(seq1, reverse_complement(seq2))
    return direct_distance <= 1 or reverse_distance <= 1


def merge_clusters(cluster1, cluster2):
    """Merges two clusters if any sequence pair between them has edit distance 0 or 1."""
    for seq1 in cluster1:
        for seq2 in cluster2:
            if is_edit_distance_one_or_zero(seq1, seq2):
                return cluster1 | cluster2
    return None


def cluster_sequences(sequences):
    """Clusters sequences based on edit distance of 0 or 1."""
    unique_sequences = list(set(sequences))
    clusters = [{seq} for seq in unique_sequences]
    while True:
        changed = False
        new_clusters = []
        merged_indices = set()

        for i, cluster1 in enumerate(clusters):
            if i in merged_indices:
                continue
            merged = False
            for j, cluster2 in enumerate(clusters):
                if i == j or j in merged_indices:
                    continue
                merged_cluster = merge_clusters(cluster1, cluster2)
                if merged_cluster:
                    new_clusters.append(merged_cluster)
                    merged_indices.update({i, j})
                    merged = True
                    changed = True
                    break
            if not merged:
                new_clusters.append(cluster1)
        clusters = new_clusters
        if not changed:
            break
    return clusters


def create_cluster_dict(clusters):
    """Creates a dictionary mapping sequences to cluster numbers."""
    sequence_to_cluster = {}
    for cluster_num, cluster in enumerate(clusters):
        for sequence in cluster:
            sequence_to_cluster[sequence] = cluster_num
    return sequence_to_cluster


def group_data_by_start(data):
    """Groups the data by the 'Start' column."""
    grouped_data = defaultdict(list)
    for row in data:
        grouped_data[row['Start']].append(row)
    return grouped_data


def remap_and_generate_outputs(data, grouped_data, fasta_output, csv_output):
    """
    Remap cluster integers to start from 1 based on encounter order, write FASTA and updated CSV files.

    Args:
        data (list): Original CSV data as a list of dictionaries.
        grouped_data (dict): Grouped data by 'Start'.
        fasta_output (str): Path to the output FASTA file.
        csv_output (str): Path to the output CSV file.
    """
    # Create a mapping table for clusters
    old_to_new = {}
    current_index = 1

    # Sequentially map cluster numbers as they are first encountered
    for start, rows in grouped_data.items():
        for row in rows:
            old_cluster = int(row['Cluster'])
            if old_cluster not in old_to_new:
                old_to_new[old_cluster] = current_index
                current_index += 1

    # Add the Spacer_Cluster_Index column to the data
    for row in data:
        row['Spacer_Cluster_Index'] = old_to_new[int(row['Cluster'])]

    # Write the updated CSV file
    with open(csv_output, mode='w', newline='') as outfile:
        fieldnames = list(data[0].keys())
        csv_writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        csv_writer.writeheader()
        csv_writer.writerows(data)

    # Write the FASTA file
    with open(fasta_output, mode='w') as fasta_file:
        for start, rows in grouped_data.items():
            # Create header
            header_parts = [
                rows[0]['Accession Number'],
                rows[0]['Start'],
                rows[0]['End'],
                rows[0]['Category']
            ]
            header = ".".join(header_parts)

            # Write header
            fasta_file.write(f">{header}\n")

            # Write remapped cluster numbers
            cluster_indices = [old_to_new[int(row['Cluster'])] for row in rows]
            fasta_file.write(", ".join(map(str, sorted(cluster_indices))) + "\n")

    print(f"FASTA file saved to {fasta_output}")
    print(f"Updated CSV file saved to {csv_output}")


def process_csv_to_fasta(input_csv, fasta_output, csv_output):
    """Main function to process the CSV, cluster sequences, and generate outputs."""
    # Step 1: Read CSV
    with open(input_csv, mode='r') as infile:
        csv_reader = csv.DictReader(infile)
        data = list(csv_reader)

    # Step 2: Extract spacer sequences
    spacer_sequences = [row['Spacer Sequence'] for row in data]

    # Step 3: Cluster the sequences
    clusters = cluster_sequences(spacer_sequences)
    sequence_to_cluster = create_cluster_dict(clusters)

    # Step 4: Add Cluster column to the data
    for row in data:
        row['Cluster'] = sequence_to_cluster[row['Spacer Sequence']]

    # Step 5: Group data by 'Start'
    grouped_data = group_data_by_start(data)

    # Step 6: Remap clusters and write FASTA/CSV outputs
    remap_and_generate_outputs(data, grouped_data, fasta_output, csv_output)


# Example Usage
if __name__ == "__main__":
    input_csv = 'original.csv'  # Input CSV file
    fasta_output = 'remapped_clusters.fasta'  # Output FASTA file
    csv_output = 'updated_clusters.csv'  # Output CSV file
    process_csv_to_fasta(input_csv, fasta_output, csv_output)
