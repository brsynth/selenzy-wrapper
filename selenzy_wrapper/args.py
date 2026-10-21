from argparse import ArgumentParser
from os import (
    path as os_path,
)

DEFAULT_NB_TARGETS = 20
DEFAULT_NB_IDS = -1
DEFAULT_HOST = "83333"
DEFAULT_MAX_NB_GENES = 5
__PACKAGE_FOLDER = os_path.dirname(os_path.realpath(__file__))
DEFAULT_DATA_FOLDER = __PACKAGE_FOLDER


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument("pathway_file", type=str, help="Path to pathway file")
    parser.add_argument("outfile", type=str, help="Path to output file")
    # parser.add_argument('--datadir',
    #                     default=None,
    #                     help='specify data directory for required databases files, please end with slash')
    parser.add_argument(
        "--nb_targets",
        type=int,
        default=DEFAULT_NB_TARGETS,
        help="Number of targets to display in results (before taxon IDs filtering) [default = 20]",
    )
    parser.add_argument(
        "--to_csv",
        type=str,
        default=None,
        help="Path to output file where genes IDs will be exported into CSV file",
    )
    parser.add_argument(
        "--nb_ids",
        type=int,
        default=DEFAULT_NB_IDS,
        help="Number of enzyme IDs to display in results (after taxon IDs filtering) [default = -1 (no limit)]",
    )
    parser.add_argument(
        "--d",
        type=float,
        default=0,
        help="Use similiarity values for preferred reaction direction only [default=0 (OFF)]",
    )
    parser.add_argument(
        "--NoMSA", action="store_true", help="Do not compute MSA/conservation scores"
    )
    parser.add_argument(
        "--smarts", action="store_true", help="Input is a reaction SMARTS string"
    )
    parser.add_argument(
        "--host_taxID",
        type=str,
        default=DEFAULT_HOST,
        help="""Taxon ID of the chassis (default: 83333 [E. coli K12])""",
    )
    parser.add_argument(
        "--enzyme_taxIDs",
        type=str,
        default=None,
        help="""Comma separated taxon IDs of output enzyme sequences""",
    )
    parser.add_argument(
        "--data-path",
        type=str,
        default=DEFAULT_DATA_FOLDER,
        help="""Path of the data directory. If does not exist,
        a 'data' folder will be created (default in package directory)""",
    )
    return parser
