from os import (
    path as os_path,
)
from argparse import (
    ArgumentParser,
    Namespace
)
from typing import (
    Dict,
)
from logging import Logger, getLogger
from json import loads as json_loads
from csv import reader as csv_reader
from pandas import read_csv as pd_read_csv
from tempfile import TemporaryDirectory
from colored import fg, bg, attr
from brs_utils import download_and_extract_tar_gz
from .Args import build_args_parser
from selenzy.Selenzy import (
    seqScore,
    updateScore,
    readData
)
from selenzy.newtax import newtax
from rptools.rplibs import (
    rpPathway
)

__SELENZY_FOLDER = 'selenzy'
__DATA_URL = 'https://gitlab.com/breakthewall/rrcache-data/-/raw/master/selenzy/data.tar.gz'
__DATA_FOLDER = os_path.join(
    __SELENZY_FOLDER,
    'data'
)

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
        '{color}{typo}rptools {version}{rst}{color} ({prog}){rst}\n'.format(
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

    if not os_path.exists(__DATA_FOLDER):
        logger.info(f'Downloading databases into {__DATA_FOLDER}...')
        download_and_extract_tar_gz(
            __DATA_URL,
            __SELENZY_FOLDER
        )
    logger.info('Reading databases...')
    pc = readData(__DATA_FOLDER)

    for rxn_id, rxn in pathway.get_reactions().items():
        try:
            with TemporaryDirectory() as tmpOutputFolder:
                # tmpOutputFolder = f'out_{rxn_id}'
                # print(rxn_id, rxn.get_smiles())
                uniprotID_score = Selenzy_score(
                    smarts=args.smarts,
                    rxn=rxn.get_smiles(),
                    host=args.host,
                    datadir=__DATA_FOLDER,
                    outdir=tmpOutputFolder,
                    nb_targets=args.nb_targets,
                    pc=pc,
                    logger=logger
                )
            # uniprotID_score_restricted = {}
            # for uniprot in uniprotID_score:
                # try:
                #     if uniprot_aaLenght[uniprot]>int(min_aa_length):
                #         uniprotID_score_restricted[uniprot] = uniprotID_score[uniprot]
                # except KeyError:
                #     logging.warning('Cannot find the following UNIPROT '+str(uniprot)+' in uniprot_aaLenght')
            rxn.add_miriam('uniprot', [i for i in uniprotID_score])
            rxn.set_selenzy_scores(uniprotID_score)
        except ValueError:
            logger.warning(f'Problem with retreiving the selenzyme information for pathway {pathway.get_id()}')

    pathway.to_rpSBML().write_to_file(args.outfile)

    logger.info(f'Results written in file \'{args.outfile}\'')


def Selenzy_score(
    smarts:  str,
    rxn: str,
    host: str,
    datadir: str,
    outdir: str,
    nb_targets: str,
    pc,
    logger: Logger = getLogger(__name__)
) -> Dict:

    uniprotID_score = {}

    newtax(
        smarts=smarts,
        smartsfile=None,
        rxn=rxn,
        host=host,
        datadir=datadir,
        outdir=outdir,
        targ=nb_targets,
        pc=pc
    )

    # ## READ AND COMPUTE SCORE
    # with open(os_path.join(
    #         outdir,
    #         'new.csv'
    #     )) as csv_file:
    #     df = pd_read_csv(csv_file, delimiter=',')
    #     for index, row in df.iterrows():
    #         uniprotID_score[row['Seq. ID']] = (
    #             100.0*row['Rxn Sim.']
    #             + 1.0*row['Consv. Score']
    #             + (-1.0)*row['Tax. distance']
    #             + (-0.1)*row['Uniprot protein evidence']
    #         )

    score = seqScore()
    data = updateScore(
        os_path.join(
            outdir,
            'new.csv'
        ),
        score
    )
    for index, row in data.iterrows():
        uniprotID_score[row['Seq. ID']] = (
            100.0*row['Rxn Sim.']
            + 1.0*row['Consv. Score']
            + (-1.0)*row['Tax. distance']
            + (-0.1)*row['Uniprot protein evidence']
    )

    # val = json_loads(data.to_json())
    # if 'Seq. ID' in val and len(val['Seq. ID'])>0:
    #     for ix in sorted(val['Seq. ID'], key=lambda z: int(z)):
    #         uniprotID_score[val['Seq. ID'][ix]] = val['Score'][ix]
    # else:
    #     raise ValueError

    return uniprotID_score


# def read_data(
#     filename: str,
#     logger: Logger = getLogger(__name__)
# ) -> Dict:
#     with open(filename, 'r') as csv_file:
#         reader = csv_reader(csv_file, delimiter=',')
#         for row in reader:
#             print(row)

if __name__ == '__main__':
    entry_point()
