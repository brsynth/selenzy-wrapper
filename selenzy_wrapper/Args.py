from argparse  import ArgumentParser
from typing import (
    Callable,
)
from selenzy_wrapper._version import __version__
from brs_utils import add_logger_args

DEFAULT_NB_TARGETS = 20
DEFAULT_NB_IDS = -1
DEFAULT_HOST = '83333'
DEFAULT_MAX_NB_GENES = 5

def build_args_parser(
    prog: str,
    description: str = '',
    epilog: str = '',
    m_add_args: Callable = None,
) -> ArgumentParser:

    parser = ArgumentParser(
        prog = prog,
        description = description,
        epilog = epilog
    )

    # Build Parser with rptools common arguments
    parser = _add_arguments(parser)

    # Add module specific arguments
    if m_add_args is None:
        parser = add_arguments(parser)
    else:
        parser = m_add_args(parser)


    return parser


def _add_arguments(parser: ArgumentParser) -> ArgumentParser:
    # Add arguments related to the logger
    parser = add_logger_args(parser)

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )

    return parser

def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument(
        'pathway_file',
        type=str,
        help='Path to pathway file'
    )
    parser.add_argument(
        'outfile',
        type=str,
        help='Path to output file'
    )
    # parser.add_argument('--datadir',
    #                     default=None,
    #                     help='specify data directory for required databases files, please end with slash')
    parser.add_argument(
        '--nb_targets',
        type=int,
        default=DEFAULT_NB_TARGETS,
        help='Number of targets to display in results (before taxon IDs filtering) [default = 20]'
    )
    parser.add_argument(
        '--to_csv',
        type=str,
        default=None,
        help='Path to output file where genes IDs will be exported into CSV file'
    )
    parser.add_argument(
        '--nb_ids',
        type=int,
        default=DEFAULT_NB_IDS,
        help='Number of enzyme IDs to display in results (after taxon IDs filtering) [default = -1 (no limit)]'
    )
    parser.add_argument(
        '--d',
        type=float,
        default=0,
        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]'
    )
    parser.add_argument(
        '--NoMSA',
        action='store_true',
        help='Do not compute MSA/conservation scores'
    )
    parser.add_argument(
        '--smarts',
        action='store_true',
        help='Input is a reaction SMARTS string'
    )
    parser.add_argument(
        '--taxonIDs',
        type=str,
        default=DEFAULT_HOST,
        help='''Comma separated taxon ids [default: 83333 (E. coli K12)].
                The first taxon ID is the one of the chassis,
                following ones are taxon IDs of output enzyme sequences'''
        )
    return parser
