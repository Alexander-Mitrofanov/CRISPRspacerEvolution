import json
import csv
import os

def get_db_path():
    """
    Get the path to the database JSON file relative to the script.

    Returns:
        str: Full path to the database JSON file.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, "identify_db.json")


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence, handling ambiguous bases."""
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N'
    }
    return ''.join(complement[base] for base in reversed(sequence))

def fetch_from_db_csv(cas_type, consensus_repeat, spacer_min_length=18, dna_alphabet=True):
    """
    Fetch spacers from the JSON database given a cas type and consensus repeat sequence.

    Args:
        cas_type (str): The cas cassette type.
        consensus_repeat (str): The consensus repeat sequence.
        spacer_min_length (int): Minimum length of spacers to include. Default is 18.
        dna_alphabet (bool): If True, ensure spacers only contain characters from "AGCTN". Default is True.

    Returns:
        list: A list of dictionaries, where each dictionary contains:
            'Spacer Sequence', 'Accession Number', 'Start', 'End', 'Category', 'Cas-Gene', 'Consensus'.
            Returns an empty list if no data is found.
    """
    json_file_path = get_db_path()

    # Load the JSON data
    try:
        with open(json_file_path, 'r') as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        return []  # Return empty structure if file not found

    # Check if the cas type exists in the database
    if cas_type not in data:
        return []  # Return empty structure if cas type is not found

    # Check if the consensus repeat exists under the cas type
    if consensus_repeat not in data[cas_type]:
        # Attempt to use reverse complement of the consensus repeat
        reverse_consensus = reverse_complement(consensus_repeat)
        if reverse_consensus not in data[cas_type]:
            return []  # Return empty structure if neither match
        consensus_repeat = reverse_consensus

    # Fetch spacers and merge them
    spacers = data[cas_type][consensus_repeat]

    # Prepare the output data
    output_data = []
    valid_dna_characters = set("AGCTN")

    for spacer_group in spacers:
        acc_number = spacer_group['acc_number']
        start = spacer_group['start']
        end = spacer_group['end']
        category = spacer_group['category']
        spacer_sequences = spacer_group['spacer_sequences'].split()  # Split spacers into a list

        for spacer in spacer_sequences:
            if len(spacer) >= spacer_min_length:
                if dna_alphabet and not set(spacer).issubset(valid_dna_characters):
                    continue  # Skip spacers with invalid characters if dna_alphabet is True
                output_data.append({
                    'Spacer Sequence': str(spacer),
                    'Accession Number': str(acc_number),
                    'Start': str(start),
                    'End': str(end),
                    'Category': str(category),
                    'Cas-Gene': str(cas_type),
                    'Consensus': str(consensus_repeat)
                })

    return output_data

