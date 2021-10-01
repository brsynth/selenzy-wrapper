from os import (
    path as os_path,
)
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger
from colored import fg, bg, attr
from .Args import build_args_parser
from rptools.rplibs import (
    rpPathway
)
from .selenzy_wrapper import selenzy_pathway


def init(
    parser: ArgumentParser,
    args: Namespace
) -> Logger:
    from brs_utils import create_logger
    from selenzy_wrapper._version import __version__

    if args.log.lower() in ['silent', 'quiet'] or args.silent:
        args.log = 'CRITICAL'

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{color}{typo}{prog} {version}{rst}{color}{rst}\n'.format(
            prog = logger.name,
            version = __version__,
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    logger.debug(args)

    return logger

def entry_point():
    parser = build_args_parser(
        prog = 'selenzy-wrapper',
        description='rpSBML wrapper for selenzy tool'
    )
    args = parser.parse_args()

    logger = init(parser, args)

    pathway = rpPathway.from_rpSBML(infile=args.pathway_file)

    selenzy_pathway(
        pathway=pathway,
        taxonIDs=args.taxonIDs,
        nb_targets=args.nb_targets,
        logger=logger
    )

    pathway.to_rpSBML().write_to_file(args.outfile)

    logger.info(f'Results written in file \'{args.outfile}\'')


if __name__ == '__main__':
    entry_point()
