
from SpacerPlacerHelpers.clustering_helpers import process_csv_to_fasta
from SpacerPlacerHelpers.clustering_helpers import fasta_to_csv
from SpacerPlacerHelpers.execute_spacer_placer import execute_spacer_placer_command

def run_spacer_placer_tool(input_fasta_file_sp, flag_use_db_sp, flag_cluster_similar_sp, result_folder_sp):
    fasta_to_csv(input_fasta_file_sp, "tmp_scv_file.csv")
    process_csv_to_fasta("tmp_scv_file.csv", "fasta_for_spacer_placer.fasta",
                         "csv_for_spacer_placer.csv")
    execute_spacer_placer_command("fasta_for_spacer_placer.fasta", result_folder_sp)


