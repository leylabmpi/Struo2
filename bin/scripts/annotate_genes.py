#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import logging
import argparse
from pprint import pprint

desc = 'Annotate a genome based on mapping to UniRef via diamond'
epi = """DESCRIPTION:
Adding UniRef IDs and taxonomy to genes from a particular genome.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('diamond_hits', metavar='diamond_hits', type=str,
                    help='tab-delim table of diamond hits')
parser.add_argument('genes_fasta_nuc', metavar='genes_fasta_nuc', type=str,
                    help='Genes in nucletide fasta format')
parser.add_argument('genes_fasta_AA', metavar='genes_fasta_AA', type=str,
                    help='Genes in amino acid fasta format')
parser.add_argument('taxonomy', metavar='taxonomy', type=str,
                    help='Taxonomy of the genome')
parser.add_argument('taxID', metavar='taxID', type=str,
                    help='NCBI TaxID of the genome')
parser.add_argument('--columns', type=str, default='qseqid,sseqid,pident,length,qstart,qend,qlen,sstart,send,slen,evalue',
                        help='Diamond output columns (default:  %(default)s)')                    
parser.add_argument('--outdir', type=str, default='genes_annotated',
                        help='Output directory (default:  %(default)s)')
parser.add_argument('--dmnd-db', type=str, default='/ebio/abt3_projects/databases_no-backup/humann2/uniref50/uniref50_annotated.1.1.dmnd',
                    help='UniRef dmnd db for annotating genes (default: %(default)s)')
parser.add_argument('--percid', type=float, default=50.0,
                        help='Percent sequence ID cutoff for calling a hit (default:  %(default)s)')
parser.add_argument('--overlap', type=float, default=80.0,
                        help='Perc. overlap cutoff (longest sequence) for calling a hit (default:  %(default)s)')
parser.add_argument('--skip', action='store_true', default=False,
                        help='Skip diamond-based annotation if the diamond hits file exists (default:  %(default)s)')
parser.add_argument('--gzip', action='store_true', default=False,
                        help='gzip output (default:  %(default)s)')
parser.add_argument('--threads', type=int, default=1,
                       help='Threads used for diamond (default:  %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def make_best_hit_index(dmnd_hit_file, outfmt_cols):
    """Making index of best diamond hits for each query
    """
    # column index
    column_idx = {x:i for i,x in enumerate(outfmt_cols)}
    # streaming hits
    hits = {}
    with open(dmnd_hit_file) as inF:
        longest_seq_len = None
        for line in inF:
            line = line.rstrip().split('\t')
            # Perc ID >= cutoff?
            try:
                pident = float(line[column_idx['pident']])
            except KeyError:
                raise KeyError('Cannot find "pident" column')
            if pident < args.percid:
                continue
            # overlap >= cutoff?
            ## longest sequence?
            try:
                query_len = float(line[column_idx['qlen']])
            except KeyError:
                raise KeyError('Cannot find "qlen" column')
            try:
                subject_len = float(line[column_idx['slen']])
            except KeyError:
                raise KeyError('Cannot find "slen" column')
            if query_len >= subject_len:
                longest_seq_len = query_len
            else:
                longest_seq_len = subject_len
            ## overlap vs longest sequence
            try:
                aln_len = float(line[column_idx['length']])
            except KeyError:
                raise KeyError('Cannot find "length" column')
            perc_overlap = aln_len / longest_seq_len * 100.0
            if perc_overlap < args.overlap:
                continue
            # better than current best hit for query?
            try:
                qseqid = line[column_idx['qseqid']]
            except KeyError:
                raise KeyError('Cannot find "qseqid" column')
            try:
                sseqid = line[column_idx['sseqid']]
            except KeyError:
                raise KeyError('Cannot find "sseqid" column')
            try:
                best_hit = hits[qseqid]
            except KeyError:
                best_hit = None
            if best_hit is None or (pident >= best_hit[1] and perc_overlap >= best_hit[2]):
                hits[qseqid] = [sseqid, pident, perc_overlap]
    return hits

def rename_seqs(best_hits, fasta_file, taxonomy, outfile, gzip_output=False):
    """Renaming sequences based on uniref hits.
    Using naming format: `gene_family|gene_length|taxonomy`
    Taxonomy format: `g__{genus};s__{species}_taxID{taxID}`
    (see https://bitbucket.org/biobakery/humann2/wiki/Home).
    """
    seq_name = None
    seq = ''
    annot_cnt = 0
    annot_skip_cnt = 0
    if gzip_output == True:
        _open = lambda x: gzip.open(x, 'ab')
        outfile += '.gz'
    else:
        _open = lambda x: open(x, 'a')
    
    with open(fasta_file) as inF, _open(outfile) as outF:
        for line in inF:
            if line.startswith('>'):
                # previous sequence
                if seq_name is not None and seq != '':
                    seq = seq.rstrip().strip('*')
                    x = '\n'.join(['>' + seq_name, seq]) + '\n'
                    if gzip_output == True:
                        x = x.encode()   
                    outF.write(x)
                    annot_cnt += 1
                else:
                    annot_skip_cnt += 1
                seq = ''
                # hit for sequence?
                query = line.rstrip().lstrip('>').split(' ')[0]
                try:
                    best_hit = best_hits[query]
                except KeyError:
                    best_hit = None
                # renaming
                if best_hit is None:
                    seq_name = None
                else:
                    seq_name = '|'.join([best_hit[0], taxonomy])
            else:
                seq += line.rstrip()
        # final sequence
        if seq_name is not None:
            seq = seq.rstrip().strip('*')
            x = '\n'.join(['>' + seq_name, seq]) + '\n'
            if gzip_output == True:
                x = x.encode()
            outF.write(x)
            annot_cnt += 1
        else:
            annot_skip_cnt += 1

    logging.info('Number of genes with an annotation: {}'.format(annot_cnt))            
    logging.info('Number of genes skipped due to no annotation: {}'.format(annot_skip_cnt))  
    logging.info('File written: {}'.format(outfile))
     
def format_taxonomy(tax, taxID):
    """Formatting taxonomy string
    """
    logging.info('Taxonomy string provided {}'.format(tax))
    logging.info('TaxID provided {}'.format(taxID))
    
    try:
        taxID = int(float(taxID.strip()))
    except ValueError:
        msg = 'ERROR: taxID "{}" is not an integer!'
        raise ValueError(msg)
    tax = re.sub('[^A-Za-z0-9-_;]+', '_', tax).split(';')
    
    if not len(tax) == 7:
        species = 's__unclassified'
    else:
        species = tax[6]
    if not len(tax) >= 6:
        genus = 'g__unclassified'
    else:
        genus = tax[5]

    if genus.startswith('G__'):
        genus = genus[3:]
    if not genus.startswith('g__'):
        genus = 'g__' + genus

    if species.startswith('S__'):
        species = genus[3:]
    if not species.startswith('s__'):
        species = 'g__' + species

    if genus == 'g__':
        genus = 'g__unclassified'
    if species == 's__':
        species = 's__unclassified'
        
    tax = '.'.join([genus, species])
    tax += '__taxID{}'.format(taxID)
    logging.info('Converted taxonomy string to {}'.format(tax))
    return tax

def main(args):
    # formatting taxonomy
    args.taxonomy = format_taxonomy(args.taxonomy, args.taxID)
    
    # filtering diamond hits
    logging.info('Finding best hit for each gene')
    outfmt_cols = args.columns.split(',')
    best_hits = make_best_hit_index(args.diamond_hits, outfmt_cols)
    #pprint(hits)

    logging.info('Renaming genes')
    # nuc
    outfile = os.path.join(args.outdir, 'annot.fna')
    rename_seqs(best_hits, args.genes_fasta_nuc,
                args.taxonomy, outfile=outfile,
                gzip_output=args.gzip)
    # AA
    outfile = os.path.join(args.outdir, 'annot.faa')
    rename_seqs(best_hits, args.genes_fasta_AA,
                args.taxonomy, outfile=outfile,
                gzip_output=args.gzip)

    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
