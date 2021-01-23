#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Uncompress gzip\'ed or bz2\'ed file'
epi = """DESCRIPTION:
Output written to STDOUT
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('input_file', metavar='input_file', type=str,
                    help='Input file')
parser.add_argument('--version', action='version', version='0.0.1')


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
    with _open(args.input_file) as inF:
        for line in inF:
            print(_decode(line, args.input_file).rstrip())
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
