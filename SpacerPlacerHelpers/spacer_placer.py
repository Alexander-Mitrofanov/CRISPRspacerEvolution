import csv
import sys

from components_spacer_placer.clustering_helpers import process_csv_to_fasta
def run_spacer_placer(input_file, flag_use_db, flag_cluster_similar):
    with open(input_file, mode='r') as file:
        csv_reader = csv.DictReader(file)

        # Parse the rows
        data = [row for row in csv_reader]

    if not flag_use_db:
        process_csv_to_fasta(input_file, "t_input/fasta_for_spacer_placer.fasta", "t_input/csv_for_spacer_placer.csv")


if __name__ == '__main__':
    input_file = '/home/alex/GitHub/CRISPRspacerEvolution/CRISPRspacerEvolution/SpacerPlacerHelpers/t_input/CAS-I-C_GTTTCAATCCACGCCCCCCGTCACCGAGGGGCGATGC.csv'
    flag_use_db = False
    run_spacer_placer(input_file, False, False )