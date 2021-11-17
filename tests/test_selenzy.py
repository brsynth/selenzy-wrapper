from unittest import TestCase
from os import path as os_path

from rptools.rplibs import rpPathway
from selenzy_wrapper.selenzy_wrapper import (
    selenzy_pathway
)


class Test_Converter(TestCase):

    folder = os_path.join(
        os_path.dirname(
            os_path.realpath(__file__)
        ),
        'data'
    )
    input = os_path.join(folder, 'lycopene.xml')
    sorted_ids = {
        'rxn_1': {
            'P21684': {
                'score': 91.109,
                'target_ID': 553
            }
        },
        'rxn_2': {
            'P21683': {
                'score': 91.9,
                'target_ID': 553
            }
        },
        'rxn_3': {
            'P21685': {
                'score': 81.82,
                'target_ID': 553
            }
        }
    }

    def test_selenzy_pathway(self):
        pathway = rpPathway.from_rpSBML(infile=self.input)
        sorted_ids = selenzy_pathway(
            pathway=pathway,
            taxonIDs='511145,553',
            nb_targets=500,
            nb_ids=1
        )
        self.assertDictEqual(
            sorted_ids,
            self.sorted_ids
        )
