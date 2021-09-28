from sys import (
    path as sys_path,
)
from os import (
    path as os_path,
)
from typing import (
    Dict
)
from logging import (
    Logger,
    getLogger
)
from tempfile import TemporaryDirectory
from pandas import read_csv as pd_read_csv
from brs_utils import download_and_extract_tar_gz
from rptools.rplibs import (
    rpPathway
)
from .Args import (
    DEFAULT_NB_TARGETS,
    DEFAULT_HOST
)
__PACKAGE_FOLDER = os_path.dirname(
    os_path.realpath(__file__)
)
__SELENZY_FOLDER = 'selenzy'
sys_path.insert(
    0,
    os_path.join(
        __PACKAGE_FOLDER,
        __SELENZY_FOLDER
    )
)
from .selenzy.Selenzy import readData
from .selenzy.newtax import newtax


__DATA_URL = 'https://gitlab.com/breakthewall/rrcache-data/-/raw/master/selenzy/data.tar.gz'
__DATA_FOLDER = os_path.join(
    __PACKAGE_FOLDER,
    'data'
)

def selenzy_pathway(
    pathway: rpPathway = None,
    taxonIDs: str = DEFAULT_HOST,
    nb_targets: int = DEFAULT_NB_TARGETS,
    logger: Logger = getLogger(__name__)
) -> None:
    if not os_path.exists(__DATA_FOLDER):
        logger.info(f'Downloading databases into {__DATA_FOLDER}...')
        download_and_extract_tar_gz(
            __DATA_URL,
            __PACKAGE_FOLDER
        )
    logger.info('Reading databases...')
    pc = readData(__DATA_FOLDER)

    for rxn_id, rxn in pathway.get_reactions().items():
        try:
            with TemporaryDirectory() as tmpOutputFolder:
                # tmpOutputFolder = f'out_{rxn_id}'
                # print(rxn_id, rxn.get_smiles())
                infos = Selenzy_infos(
                    smarts=True,
                    rxn=rxn.get_smiles(),
                    taxonIDs=taxonIDs,
                    datadir=__DATA_FOLDER,
                    outdir=tmpOutputFolder,
                    nb_targets=nb_targets,
                    pc=pc,
                    logger=logger
                )
            # Sort descending order
            sorted_infos = dict(
                sorted(
                    infos.items(),
                    key=lambda item: item[1]['score'],
                    reverse=True
                )
            )

            # uniprotID_score_restricted = {}
            # for uniprot in uniprotID_score:
                # try:
                #     if uniprot_aaLenght[uniprot]>int(min_aa_length):
                #         uniprotID_score_restricted[uniprot] = uniprotID_score[uniprot]
                # except KeyError:
                #     logging.warning('Cannot find the following UNIPROT '+str(uniprot)+' in uniprot_aaLenght')
            rxn.add_miriam('uniprot', [i for i in sorted_infos])
            rxn.set_selenzy_infos(sorted_infos)
        except ValueError:
            logger.warning(f'Problem with retreiving the selenzyme information for pathway {pathway.get_id()}')


def Selenzy_infos(
    smarts:  str,
    rxn: str,
    taxonIDs: str,
    datadir: str,
    outdir: str,
    nb_targets: str,
    pc,
    logger: Logger = getLogger(__name__)
) -> Dict:

    infos = {}

    newtax(
        smarts=smarts,
        smartsfile=None,
        rxn=rxn,
        host=taxonIDs,
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

    # score = seqScore()
    # data = updateScore(
    #     os_path.join(
    #         outdir,
    #         'new.csv'
    #     ),
    #     score
    # )

    # Split taxon IDs
    taxonIDs = taxonIDs.split(',')

    with open(os_path.join(
            outdir,
            'new.csv'
        )) as csv_file:
        df = pd_read_csv(csv_file, delimiter=',')
        # If a single taxon ID is passed,
        # then return everything
        if len(taxonIDs) == 1:
            for index, row in df.iterrows():
                infos[row['Seq. ID']] = set_infos(
                    score=row['Score'],
                    target_id=row['target_ID'],
                    logger=logger
                )
        # Otherwise, return only those have
        # negative taxonomic distance 
        else:
            for index, row in df.iterrows():
                # Ignore the host taxon ID
                for taxonID in taxonIDs[1:]:
                    # If a negative distance is found,
                    # then store the seqID and
                    # pass to the next line.
                    # Otherwise, check next taxonID distance
                    if row[f'{taxonID}_target_host_dist'] <= 0:
                        infos[row['Seq. ID']] = set_infos(
                            score=row['Score'],
                            target_id=row['target_ID'],
                            logger=logger
                        )
                        break

    # val = json_loads(data.to_json())
    # if 'Seq. ID' in val and len(val['Seq. ID'])>0:
    #     for ix in sorted(val['Seq. ID'], key=lambda z: int(z)):
    #         uniprotID_score[val['Seq. ID'][ix]] = val['Score'][ix]
    # else:
    #     raise ValueError

    return infos


def set_infos(
    score: float,
    target_id: str,
    logger: Logger = getLogger(__name__)
) -> Dict:
    return  {
        'score': score,
        'target_ID': target_id
    }

