# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
from os.path import join
from shutil import rmtree


in_file_pathway = Path(__file__).resolve().parent / 'data' / 'input' / 'pathway.xml'


def test_output(tmpdir):

    #compare the contents of output file actual vs expected
    #python -m selenzy_wrapper 'tests/data/input/pathway.xml' 'tests/data/output/output_test.xml' --nb_targets '500' --d '0.0' --taxonIDs '83333' --to_csv 'tests/data/output/output_test.csv'

    result = subprocess.run(
        [
            'python','-m', 'selenzy_wrapper',
            in_file_pathway,
            join(tmpdir, "output_test.xml"),
            '--nb_targets', '500',
            '--d', '0.0',
            '--taxonIDs', '83333',
            '--to_csv', join(tmpdir, "output_test.csv")
        ],
        stdout=subprocess.PIPE
    )

    out_file = tmpdir.join("output_test.xml").strpath
    with open(out_file) as f:
        actual = f.readlines()
        actual.sort()

    fname = Path(__file__).resolve().parent / 'data' / 'output' / 'output_test.xml'
    with open(fname) as f:
        expected = f.readlines()
        expected.sort()

    print('actual:', actual)
    print('expected:', expected)
    assert actual == expected

    out_file = tmpdir.join("output_test.csv").strpath
    with open(out_file) as f:
        actual = f.readlines()
        actual.sort()

    fname = Path(__file__).resolve().parent / 'data' / 'output' / 'output_test.csv'
    with open(fname) as f:
        expected = f.readlines()
        expected.sort()

    print('actual:', actual)
    print('expected:', expected)
    assert actual == expected

    rmtree(tmpdir)
