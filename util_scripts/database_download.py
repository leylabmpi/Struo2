#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
import functools
import multiprocessing as mp
# 3rd party
import requests
import bs4

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Download Struo2 database files'
epi = """DESCRIPTION:
A helper script for downloading pre-built custom database files.
Multiple GTDB releases & databases (eg., kraken2 or humann3) 
can be downloaded.

Note: use "--" to separate --databass parameters from <output_dir>
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('output_dir', type=str,
                    help='Output directory')
parser.add_argument('-r', '--release', type=str, nargs='+',
                    choices = ['95', '202'], default=['202'],
                    help='GTDB release')
parser.add_argument('-d', '--database', type=str, nargs='+',
                    choices = ['kraken2', 'humann3', 'taxdump', 'phylogeny',
                               'metadata', 'genes'], default=['metadata'],
                    help='Database(s) to download ')
parser.add_argument('-u', '--base-url', type=str,
                    default='http://ftp.tue.mpg.de/ebio/projects/struo2/',
                    help='Base url for downloads')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Parallel download processes')
parser.add_argument('--version', action='version', version='0.0.1')

# functions
def decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def dl_file(l, url, out_dir):
    """
    Download file from url
    """
    if l['href'].startswith('?') or l['href'].endswith('/'):
        return None
    with requests.get(url + '/' + l['href'].lstrip('/'), stream=True) as r:
        if r.status_code != 404:
            comp = os.path.splitext(l['href'])[1] in ('.gz', '.bz2')
            out_file = os.path.join(out_dir, l['href'])
            if comp:
                with open(out_file, 'wb') as outF:
                    for chunk in r.iter_content(chunk_size = 1024):
                        if chunk:
                            outF.write(chunk)
            else:
                with open(out_file, 'w') as outF:
                    for line in r.iter_lines(decode_unicode=True):
                        outF.write(decode(line) + '\n')                        
            logging.info('File written: {}'.format(out_file))
    return None
            
def dl_files(base_url, release, database, out_dir, threads):
    """
    List files from url and download all available
    """
    # output directory
    out_dir = os.path.join(out_dir, release, database)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    # base url: GET
    url = os.path.join(base_url, release, database)
    r = requests.get(url)
    if r.status_code == 404:
        logging.warning('WARNING: 404 status code for url: {}'.format(url))
        return None
    # file urls: GET
    data = bs4.BeautifulSoup(r.text, 'html.parser')
    func = functools.partial(dl_file, url=url, out_dir=out_dir)
    if args.threads > 1:
        pool = mp.Pool(threads)
        pool.map(func, data.find_all('a'))
    else:
        [x for x in map(func, data.find_all('a'))]
    try:
        pool.close()
    except UnboundLocalError:
        pass
        
def main(args):
    # args
    args.release = ['GTDB_release{}'.format(x) for x in args.release]
    # output
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    # list files
    for db in args.database:
        for release in args.release:
            dl_files(args.base_url, release, db,
                     args.output_dir, args.threads)
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
