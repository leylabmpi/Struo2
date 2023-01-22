#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import resource
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

Note: use "--" to separate "--database" parameters from <output_dir>
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('output_dir', type=str,
                    help='Output directory')
parser.add_argument('-r', '--release', type=str, nargs='+',
                    choices = ['95', '202', '207'], default=['207'],
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
parser.add_argument('-m', '--max-recursion', type=int, default=1048576,
                    help='Max recursion limit')
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

def write_lines(url, l, out_dir):
    """"
    Writing lines obtained from requests 
    """
    with requests.get(url + '/' + l['href'].lstrip('/'), stream=True) as r:        
        if r.status_code == 404:
            return None
        if l['href'] == 'database.kraken':  # debug
            return None
        out_file = os.path.join(out_dir, l['href'])
        with open(out_file, 'w') as outF:
            for line in r.iter_lines(decode_unicode=True):
                outF.write(decode(line) + '\n')
    logging.info(f'File written: {out_file}')
    
def write_chunks(url, l, out_dir):
    """
    Writing chunks obtained from requests
    """
    with requests.get(url + '/' + l['href'].lstrip('/'), stream=True) as r:        
        if r.status_code == 404:
            return None
        out_file = os.path.join(out_dir, l['href'])
        with open(out_file, 'wb') as outF:
            for chunk in r.iter_content(chunk_size = 1024):
                if chunk:
                    outF.write(chunk)
    logging.info(f'File written: {out_file}')

def dl_file(l, url, out_dir):
    """
    Download file from url
    """
    if l['href'].startswith('?') or l['href'].endswith('/'):
        return None
    try:
        write_lines(url, l, out_dir)
    except UnicodeDecodeError:
        write_chunks(url, l, out_dir)    
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
    bs4_list = [x for x in data.find_all('a')]
    if args.threads > 1:
        pool = mp.Pool(threads)
        pool.map(func, bs4_list)
    else:
        [x for x in map(func, bs4_list)]
    try:
        pool.close()
    except UnboundLocalError:
        pass
        
def set_recursion(max_rec):
    """
    max_rec = 0x100000
    """
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)
    logging.info(f'Max recursion set to: {max_rec}')
    
def main(args):
    # args
    args.release = ['GTDB_release{}'.format(x) for x in args.release]
    # recursion
    set_recursion(args.max_recursion)
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
