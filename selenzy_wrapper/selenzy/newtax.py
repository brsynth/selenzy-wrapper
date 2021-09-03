""" 
usage: newtax.py [-h] [-tar TAR] [-d D] [-outfile OUTFILE] [-NoMSA] [-smarts]
                 [-smartsfile] [-host HOST]
                 rxn datadir outdir

Run Selenzyme for multiple hosts 

positional arguments:
  rxn               Input reaction [default = rxn file]
  datadir           specify data directory for required databases files,
                    please end with slash
  outdir            specify output directory for all output files, including
                    final CSV file, please end with slash

optional arguments:
  -h, --help        show this help message and exit
  -tar TAR          Number of targets to display in results [default = 20]
  -d D              Use similiarity values for preferred reaction direction
                    only [default=0 (OFF)]
  -outfile OUTFILE  specify non-default name for CSV file output
  -NoMSA            Do not compute MSA/conservation scores
  -smarts           Input is a reaction SMARTS string
  -smartsfile       Input is a reaction SMARTS file
  -host HOST        Comma separated taxon ids [default: E. coli]
  
Example: 
    python newtax.py "O=C([O-])CCC(=O)C(=O)[O-].NC(CC(=O)[O-])C(=O)O>>O=C([O-])CC(=O)C(=O)[O-].NC(CCC(=O)[O-])C(=O)O" datadir outdir -host "83333,90370,1236,590" 

Notes:
    
This example code presents two possible modes of operation of the Selenzyme tool. 

The call to Selenzy2.analyse2 accepts both as single string host or an array with different hosts. 
In the command line version, if different hosts are to be entered, the host IDs must be separated by commas and not spaced.

A csv file with the name "result.csv" is generated, as it was generated from the original Selenzyme code, using the first taxon_ID or the default host if is a taxonomic group, i.e. Gammaproteobacteria.

If the host is a set of hosts, the code generates "new.csv", which adds the following columns to the original csv (for each host):
- Host/Target distance
- Host/Common ancestor distance
- Target/Common ancestor distance
In addition, a final column is generated where the Target ID is indicated.
"""
import csv
import pandas as pd
import Selenzy
import Selenzy2
import os
import argparse
from typing import (
    List,
)

fileSeqOrg = "seq_org.tsv" 
fileLineage = "org_lineage.csv"

def newtax(
    smarts: str,
    smartsfile: str,
    rxn: str,
    host: str,
    datadir: str,
    outdir: str,
    targ: str,
    pc = None,
):
    newpath = os.path.join(outdir)
    if not os.path.exists(newpath):
        os.makedirs(newpath) 
    if smarts is not None:
        rxnInput = ['-smarts', rxn]
    elif smartsfile:
        rxnInput = ['-smartsfile', rxn]
    else:
        rxnInput = ['-rxn', rxn]
    # Convert host list to array
    host = split_host(host)
    host_analyse(
        rxnInput=rxnInput,
        host=host,
        targ=targ,
        datadir=datadir,
        outdir=outdir,
        pc=pc
    )

def split_host(host: str) -> List[str]:
    return host.split(",")

def host_analyse(
    rxnInput: str,
    datadir: str,
    outdir: str,
    host: str,
    targ: str,
    pc = None
):
    tax = Selenzy.readTaxonomy(datadir, fileLineage)
    taxNodes = Selenzy2.superTax2(tax)
    if host == str():
        datos = Selenzy.analyse(
            rxnInput=rxnInput,
            targ=targ,
            datadir=datadir,
            outdir=outdir,
            csvfilename="result.csv",
            NoMSA=True,
            host=host
        )
    else:
        _host = host[0]
        if _host not in tax:
            _host = '83333'
        datos = Selenzy.analyse(
            rxnInput=rxnInput,
            targ=targ,
            datadir=datadir,
            outdir=outdir,
            csvfilename="result.csv",
            NoMSA=True,
            host=_host,
            pc=pc
        )
        data = Selenzy.updateScore(os.path.join(outdir,"result.csv"), Selenzy.seqScore())
        outfile = os.path.join(outdir,'new.csv')
        df2 = Selenzy2.analyse2(
            df=data,
            host=host,
            taxNodes=taxNodes,
            datadir=datadir,
            outfile=outfile
        )

def arguments():
    parser = argparse.ArgumentParser(description="""Run Selenzyme for multiple hosts.""")
    parser.add_argument(
        'rxn', 
        help='Input reaction [default = rxn file]'
    )
    parser.add_argument(
        '-tar',
        type=float,
        default=20,
        help='Number of targets to display in results [default = 20]'
        )
    parser.add_argument(
        '-d',
        type=float,
        default=0,
        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]'
    )
    parser.add_argument(
        'datadir',
        help='specify data directory for required databases files, please end with slash'
    )
    parser.add_argument(
        'outdir',
        help='specify output directory for all output files, including final CSV file, please end with slash'
    )
    parser.add_argument(
        '-outfile',
        help='specify non-default name for CSV file output'
    )
    parser.add_argument(
        '-NoMSA',
        action='store_true',
        help='Do not compute MSA/conservation scores'
    )
    parser.add_argument(
        '-smarts',
        action='store_true',
        help='Input is a reaction SMARTS string'
    )
    parser.add_argument(
        '-smartsfile',
        action='store_true',
        help='Input is a reaction SMARTS file'
    )
    parser.add_argument(
        '-host',
        type=str,
        default='83333',
        help='Comma separated taxon ids [default: E. coli]'
    )
    arg = parser.parse_args()
    return arg

if __name__ == '__main__':
    args = arguments()
    newtax(
        smarts=args.smarts,
        smartsfile=args.smartsfile,
        rxn=args.rxn,
        host=args.host,
        datadir=args.datadir,
        outdir=args.outdir,
        targ=args.tar
    )
