import json
import csv

def fetch_from_db(cas_type, consensus_repeat, json_file_path, spacer_min_length=18, dna_alphabet=True):
    """
    Fetch spacers from the JSON database given a cas type and consensus repeat sequence.

    Args:
        cas_type (str): The cas cassette type.
        consensus_repeat (str): The consensus repeat sequence.
        json_file_path (str): Path to the JSON file.
        spacer_min_length (int): Minimum length of spacers to include. Default is 18.
        dna_alphabet (bool): If True, ensure spacers only contain characters from "AGCT". Default is True.

    Returns:
        list: A list of lists, where each sublist contains:
            [index, spacer_sequence, accession_number, start, end, category]
    """
    # Load the JSON data
    try:
        with open(json_file_path, 'r') as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file at {json_file_path} was not found.")

    # Check if the cas type exists in the database
    if cas_type not in data:
        raise ValueError(f"The cas type '{cas_type}' does not exist in the database.")

    # Check if the consensus repeat exists under the cas type
    if consensus_repeat not in data[cas_type]:
        raise ValueError(f"The consensus repeat '{consensus_repeat}' does not exist for the cas type '{cas_type}'.")

    # Fetch spacers and merge them
    spacers = data[cas_type][consensus_repeat]

    # Prepare the output data
    output_data = []
    index = 1
    valid_dna_characters = set("AGCT")

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
                output_data.append([
                    index, spacer, acc_number, start, end, category
                ])
                index += 1

    return output_data


def check():
    data = fetch_from_db("CAS-I-E",
                         "GTGTTCCCCGCGCCAGCGGGGATAAACCG",
                         "identify_db.json")
    for row in data:
        print(row)


check()