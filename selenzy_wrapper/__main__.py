from brs_utils import build_args_parser
from brs_utils import init as init_logger
from pandas import DataFrame
from rplibs import rpPathway

from ._version import __version__
from .args import DEFAULT_NB_IDS, add_arguments
from .selenzy_wrapper import selenzy_pathway


def entry_point():
    parser = build_args_parser(
        prog="selenzy-wrapper",
        version=__version__,
        description="rpSBML wrapper for selenzy tool",
        m_add_args=add_arguments,
    )
    args = parser.parse_args()

    logger = init_logger(parser, args, __version__)

    pathway = rpPathway(args.pathway_file)

    taxonIDs = args.host_taxID
    if args.enzyme_taxIDs is not None:
        taxonIDs += f",{args.enzyme_taxIDs}"

    selenzy_pathway(
        pathway=pathway,
        taxonIDs=taxonIDs,
        nb_targets=args.nb_targets,
        nb_ids=args.nb_ids,
        datadir=args.data_path,
        logger=logger,
    )

    pathway.to_rpSBML().write_to_file(args.outfile)

    if args.to_csv is not None:
        genes = selenzinfo2table(pathway=pathway, maxgenes=args.nb_ids)
        genes.to_csv(args.to_csv, index=False)

    logger.info(f"Results written in file '{args.outfile}'")


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
        rxn.get_id(): rxn.get_selenzy() for rxn in pathway.get_list_of_reactions()
    }
    genes = DataFrame(columns=["Name", "Type", "Part", "Step"])
    for rxn_id, enz_ids in selenzyme_info.items():
        nb_genes = 0
        if maxgenes == DEFAULT_NB_IDS or nb_genes < maxgenes:
            nb_genes += 1
            for enz_id in enz_ids:
                genes.loc[len(genes)] = [
                    enz_id,
                    "gene",
                    enz_id,
                    pathway.get_reaction(rxn_id).get_idx_in_path(),
                ]
    return genes


if __name__ == "__main__":
    entry_point()
