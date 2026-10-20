import csv
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree

from selenzy_wrapper.args import DEFAULT_DATA_FOLDER

DATA_FOLDER = Path(__file__).resolve().parent / "data"


def read_enzyme_annotations(path):
    namespaces = {
        "sbml": "http://www.sbml.org/sbml/level3/version1/core",
        "brsynth": "http://brsynth.eu",
    }
    root = ElementTree.parse(path).getroot()
    return {
        reaction.attrib["id"]: {
            annotation.tag: {field.tag: field.attrib for field in annotation}
            for annotation in reaction.findall(".//brsynth:brsynth/*", namespaces)
            if annotation.tag.startswith("{http://brsynth.eu}selenzy_")
        }
        for reaction in root.findall(".//sbml:reaction", namespaces)
    }


def test_output(tmp_path):
    output_xml = tmp_path / "output_test.xml"
    output_csv = tmp_path / "output_test.csv"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "selenzy_wrapper",
            str(DATA_FOLDER / "input" / "pathway.xml"),
            str(output_xml),
            "--nb_targets",
            "500",
            "--d",
            "0.0",
            "--host_taxID",
            "83333",
            "--to_csv",
            str(output_csv),
            "--data-path",
            DEFAULT_DATA_FOLDER,
        ],
        stdout=subprocess.PIPE,
        check=True,
    )

    # Compare enzyme IDs, scores and taxon IDs independently of SBML serialization.
    expected_annotations = read_enzyme_annotations(
        DATA_FOLDER / "output" / "output_test.xml"
    )
    assert expected_annotations
    assert all(expected_annotations.values())
    assert read_enzyme_annotations(output_xml) == expected_annotations

    with output_csv.open(newline="") as actual_file:
        actual = list(csv.reader(actual_file))
    with (DATA_FOLDER / "output" / "output_test.csv").open(newline="") as expected_file:
        expected = list(csv.reader(expected_file))
    assert actual[0] == expected[0]
    assert sorted(actual[1:]) == sorted(expected[1:])
