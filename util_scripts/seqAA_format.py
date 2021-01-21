#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging

desc = 'Converting AA sequence file to format used for Struo2'
epi = """DESCRIPTION:



OUTPUT
* fasta of renamed AA sequences 
* tab-delim metadata table for each AA sequence
  * columns
    * seq_uuid => UUID for sequence
    * seq_orig_name => original sequence name
    * genus => genus of the gene
    * species => species of the gene
    * taxid => taxid (species-level) of the gene
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('input_file', metavar='input_file', type=str,
                    help='Input file')
parser.add_argument('--version', action='version', version='0.0.1')
#For default: (default: %(default)s)


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
# logging.info()


def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(line, infile):
    """
    Decoding input, depending on the file extension
    """
    if os.path.isfile(infile) and (infile.endswith('.gz') or infile.endswith('.bz2')):
        line = line.decode('utf-8')
    return line



def main(args):
    pass
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
