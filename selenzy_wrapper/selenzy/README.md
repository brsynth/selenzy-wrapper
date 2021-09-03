Enzyme selection tool

usage: newtax.py [-h] [-tar TAR] [-d D] [-outfile OUTFILE] [-NoMSA] [-smarts]
                 [-smartsfile] [-host HOST] rxn datadir outdir

Run Selenzyme for multiple hosts 

Positional arguments:

  * rxn               Input reaction [default = rxn file]

  * datadir           specify data directory for required databases files,
                    please end with slash

  * outdir            specify output directory for all output files, including
                    final CSV file, please end with slash

Optional arguments: 

  * -h, --help        show this help message and exit

  * -tar TAR          Number of targets to display in results [default = 20]

  * -d D              Use similiarity values for preferred reaction direction
                    only [default=0 (OFF)]

  * -outfile OUTFILE  specify non-default name for CSV file output

  * -NoMSA            Do not compute MSA/conservation scores

  * -smarts           Input is a reaction SMARTS string

  * -smartsfile       Input is a reaction SMARTS file

  * -host HOST        Comma separated taxon ids [default: E. coli]
  
Example: 
```(bash)
    python newtax.py "O=C([O-])CCC(=O)C(=O)[O-].NC(CC(=O)[O-])C(=O)O>>O=C([O-])CC(=O)C(=O)[O-].NC(CCC(=O)[O-])C(=O)O" datadir outdir -host "83333,90370,1236,590"
```

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

