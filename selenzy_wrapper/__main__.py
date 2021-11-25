from os import (
    path as os_path,
)
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger
from colored import fg, bg, attr
from pandas import DataFrame
from .Args import (
    build_args_parser,
    DEFAULT_NB_IDS
)
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
        nb_ids=args.nb_ids,
        logger=logger
    )

    pathway.to_rpSBML().write_to_file(args.outfile)

    if args.to_csv is not None:
        genes = selenzinfo2table(
            pathway=pathway,
            maxgenes=args.nb_ids
        )
        genes.to_csv(args.to_csv, index=False)

    logger.info(f'Results written in file \'{args.outfile}\'')

def selenzinfo2table(pathway, maxgenes=DEFAULT_NB_IDS):
    """Convert the selenzyme_info dictionary into the input table: Name, Type, Part, Step

    It assumes that pathway steps are in reverse direction and are called as RP1, RP2, etc.

    :param si: The selenzyme information of the heterologous pathway
    :param maxgenes: The maximal number of genes

    :type si: dict
    :type maxgenes: int

    :rtype: pandas.DataFrame
    :return: The table of parts
    """
    selenzyme_info = {
        rxn.get_id(): rxn.get_selenzy()
        for rxn in
        pathway.get_list_of_reactions()
    }
    genes = DataFrame(columns=['Name','Type', 'Part', 'Step'])
    for rxn_id in selenzyme_info:
        nb_genes = 0
        if maxgenes == DEFAULT_NB_IDS or nb_genes < maxgenes:
            nb_genes += 1
            for enz_id in selenzyme_info[rxn_id]:
                genes.loc[len(genes)] = [
                    enz_id,
                    'gene',
                    enz_id,
                    pathway.get_reaction(rxn_id).get_idx_in_path()
                ]
    return genes


if __name__ == '__main__':
    entry_point()
