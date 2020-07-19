#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("curl"))
suppressPackageStartupMessages(library("data.table"))

# create parser object
parser <- ArgumentParser()

# specifying options
parser$add_argument("metadata_urls", nargs='+', help=">=1 url to GTDB metadata")
parser$add_argument("-o", "--output", type='character', default='metadata.tsv',
			   help="Output file name [default: %(default)s]")
parser$add_argument("-c", "--columns", type='character', default='ncbi_organism_name,ncbi_genbank_assembly_accession,scaffold_count,contig_count,gc_percentage,genome_size,checkm_completeness,checkm_contamination,checkm_strain_heterogeneity,ncbi_assembly_level,ncbi_refseq_category,ncbi_species_taxid,ncbi_taxonomy,gtdb_taxonomy,mimag_high_quality,gtdb_representative',
			   help="Table columns to keep [default: %(default)s]")
parser$add_argument("-f", "--filter", type='character', default='gtdb_representative == "t" & checkm_completeness >= 50 & checkm_contamination < 5',
			   help="Table columns to keep [default: %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
			   help="Print extra output [default: %(default)s]")
parser$add_argument("-q", "--quietly", action="store_false",
			   dest="verbose", help="Print little output")
args <- parser$parse_args()


# reading in table(s)
write(sprintf('Keeping columns: %s', args['columns']), stderr())
cols = unlist(strsplit(unlist(args['columns']), ','))
write('----', stderr())

df = list()
for(url in unlist(args['metadata_urls'])){
    write(sprintf('Reading in file: %s', url), stderr())
    df[[url]] = fread(url, sep='\t', check.names=TRUE)[, ..cols]
}

df = do.call(rbind, df)
x = as.character(nrow(df))
write(sprintf('Number of rows in the combined table: %s', x), stderr())

# Filtering
x = unlist(args['filter'])[1]
write(sprintf('Filtering rows by expression: %s', x), stderr())
df = df[eval(parse(text=x)),]
x = as.character(nrow(df))
write(sprintf('Number of rows after filtering: %s', x), stderr())

# Writing table
write(sprintf('Writing file to: %s', args['output']), stderr())
out_file = unlist(args['output'])
fwrite(df, file=out_file, sep='\t', quote=FALSE, row.names=FALSE)