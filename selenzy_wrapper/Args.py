from argparse  import ArgumentParser
from typing import (
    Callable,
)
from selenzy_wrapper._version import __version__
from typing import(
    Callable,
)
from brs_utils import add_logger_args


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
    parser.add_argument('pathway_file',
                        help='path to pathway file')
    parser.add_argument('outfile',
                        help='specify output file')
    # parser.add_argument('--datadir',
    #                     default=None,
    #                     help='specify data directory for required databases files, please end with slash')
    parser.add_argument('--nb_targets', type=int, default=20,
                        help='Number of targets to display in results [default = 20]')
    parser.add_argument('--d', type=float, default=0,
                        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]')
    parser.add_argument('--NoMSA', action='store_true',
                        help='Do not compute MSA/conservation scores')
    parser.add_argument('--smarts', action='store_true',
                        help='Input is a reaction SMARTS string')
    parser.add_argument('--host', type=str, default='83333',
                        help='Comma separated taxon ids [default: E. coli]')
    return parser
