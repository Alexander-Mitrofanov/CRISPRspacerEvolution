import pathlib
import os
import argparse

from SpacerPlacerHelpers.spacer_placer import run_spacer_placer_tool


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta', dest='fasta_file',
                        help='Fasta file path (it can be either protein or DNA, see -st and -sc for details).',
                        default='Example/NC_006513.fa')

    parser.add_argument('--tool', type=str, default="identify",
                        help='what tool is used')

    parser.add_argument('--model', type=str, default="ALL",
                       help='model_to_use (default: ALL)')

    parser.add_argument('--additional_model', type=str, default=None,
                        help='model_to_use (default: None)')

    parser.add_argument('--result_folder', type=str, default="Results",
                        help='folder with the result (default: Results)')

    parser.add_argument('--pickle_report', type=str, default='',
                        help='pickled report file (default: None)')

    parser.add_argument('--strand', type=str, default=True,
                        help='CRISPR array orientation prediction (default: True)')

    parser.add_argument('--cas', type=str, default=False,
                        help='cas genes computation (default: False)')

    parser.add_argument('--is_element', type=str, default=True,
                        help='is element computation (default: True)')

    parser.add_argument('--parallel', type=str, default=True,
                        help='parallel computations (default: True)')

    parser.add_argument('--fast_run', type=str, default=False,
                        help='fast run option (default: False)')

    parser.add_argument('--degenerated', type=bool, default=True,
                        help='degenerated_repeat_computation (default: True)')

    parser.add_argument('--min_len_rep', type=int, default=21,
                        help='min avg. length of the repeats (default: 21)')

    parser.add_argument('--max_len_rep', type=int, default=55,
                        help='max avg. length of the repeats (default: 55)')

    parser.add_argument('--min_len_spacer', type=int, default=18,
                        help='min avg. length of spacers (default: 18)')

    parser.add_argument('--max_len_spacer', type=int, default=78,
                        help='max avg. length of spacers (default: 78)')

    parser.add_argument('--min_repeats', type=int, default=3,
                        help='min number of repeats (default: 3)')

    parser.add_argument('--enhancement_max_min', type=bool, default=True,
                        help='enhancement with filter (default: True)')

    parser.add_argument('--enhancement_start_end', type=bool, default=True,
                        help='enhancement with start end omitting (default: True)')

    parser.add_argument('--max_identical_spacers', type=int, default=4,
                        help='maximum number of identical spacers in the array (default: 4)')

    parser.add_argument('--max_identical_cluster_spacers', type=int, default=3,
                        help='maximum number of consecutive identical spacers in the array (default: 3)')

    parser.add_argument('--margin_degenerated', type=int, default=30,
                        help='maximum length of the spacer margin for the degenerated search (default: 30)', )

    parser.add_argument('--max_edit_distance_enhanced', type=int, default=6,
                        help='maximum edit distance for the evaluated array enhancement (default: 6)')

    parser.add_argument(
        '--input_fasta_file_sp',
        required=True,
        type=str,
        help="Path to the input FASTA file."
    )

    parser.add_argument(
        '--flag_use_db_sp',
        type=bool,
        help="Flag to use the database for Spacer Placer."
    )

    parser.add_argument(
        '--flag_cluster_similar_sp',
        type=bool,
        help="Flag to cluster similar spacer sequences."
    )

    parser.add_argument(
        '--folder_output_sp',
        required=True,
        help="Path to the folder where the output will be stored."
    )

    args = parser.parse_args()
    return args


def run_crispr_identify(args):
    """
    This function runs the CRISPRidentify tool.

    Parameters:
        The tool's parameters can be seen from the line parser = argparse.ArgumentParser()

    Returns:
        dirname_identify: Returns the result paths to future functions.
    """

    cur_path = str(pathlib.Path().absolute())
    dirname_identify = cur_path + '/tmp/output-CRISPRidentify'

    result = os.system('python3.7 CRISPRidentify/CRISPRidentify.py --file ' + args.fasta_file +
              ' --model ' + args.model +
              ' --result_folder ' + dirname_identify +
              ' --strand ' + str(args.strand) +
              ' --cas ' + "True" +
              ' --is_element ' + str(args.is_element) +
              ' --parallel ' + str(args.parallel) +
              ' --fast_run ' + str(args.fast_run) + ' --degenerated ' + str(args.degenerated) +
              ' --min_len_rep ' + str(args.min_len_rep) +
              ' --max_len_rep ' + str(args.max_len_rep) +
              ' --min_len_spacer ' + str(args.min_len_spacer) +
              ' --max_len_spacer ' + str(args.max_len_spacer) +
              ' --min_repeats ' + str(args.min_repeats) +
              ' --enhancement_max_min ' + str(args.enhancement_max_min) +
              ' --enhancement_start_end ' + str(args.enhancement_start_end) +
              ' --max_identical_spacers ' + str(args.max_identical_spacers) +
              ' --max_identical_cluster_spacers ' + str(args.max_identical_cluster_spacers) +
              ' --margin_degenerated ' + str(args.margin_degenerated) +
              ' --max_edit_distance_enhanced ' + str(args.max_edit_distance_enhanced))
    os.chdir(cur_path)
    return dirname_identify


def run_spacer_placer(args):
    input_fasta_file_sp = args.input_fasta_file_sp
    flag_use_db_sp = args.flag_use_db_sp
    flag_cluster_similar_sp = args.flag_cluster_similar_sp
    folder_output_sp = args.folder_output_sp
    run_spacer_placer_tool(input_fasta_file_sp,
                           flag_use_db_sp,
                           flag_cluster_similar_sp,
                           folder_output_sp)

def main():
    args = parse_arguments()
    if args.tool == "identify":
        run_crispr_identify(args)
    else:
        print("running sp")
        run_spacer_placer(args)


if __name__ == '__main__':
    print("here")
    main()
