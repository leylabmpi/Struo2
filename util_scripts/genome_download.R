#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

# create parser object
parser <- ArgumentParser()

# specifying options
parser$add_argument("acc_table", nargs=1, help="Table containing assembly accessions (tab-delim with header)")
parser$add_argument("-c", "--column", type='character', default='ncbi_genbank_assembly_accession',
                    help="Column name containing accessions [default: %(default)s]")
parser$add_argument("-o", "--output", type='character', default='.',
                    help="Path for output [default: %(default)s]")
parser$add_argument("-p", "--procs", type='integer', default=1,
	            help="Number of parallel processes [default: %(default)s]")
parser$add_argument("-r", "--retries", type='integer', default=3,
                    help="Number of retries [default: %(default)s]")
parser$add_argument("-d", "--database", type='character', default='genbank',
                    help="database to download (-s flag for ncbi-genome-download) [default: %(default)s]")
parser$add_argument("-x", "--params", type='character', default='archaea,bacteria',
	            help="Filtering parameters for ncbi-genome-download [default: %(default)s]")
parser$add_argument("-f", "--filter", action="store_true", default=FALSE,
                    help="Check for 'fasta_file_path' and just download any accessions lacking values [default: %(default)s]")
parser$add_argument("-s", "--skip", action="store_true", default=FALSE,
                    help="Skip the genome downloading; useful if re-running to re-make the output table [default: %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default: %(default)s]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
args = parser$parse_args()


# checking for executables
exe = 'ncbi-genome-download'
hits = unlist(Sys.which('ncbi-genome-download'))
if(hits[1] == ''){
    stop(sprintf('Cannot find executable: %s', exe))
}

# reading in table
x = unlist(args['acc_table'])[1]
write(sprintf('Reading table: %s', x), stderr())
df = read.delim(x, sep='\t')
write(sprintf('Number of rows: %s', nrow(df)), stderr())

# filtering table
## checking for "fasta_file_path" column
filter_bool = unlist(args['filter'])[1]
if(filter_bool == TRUE){
    write('Filtering to just rows with NAs in `fasta_file_path` column', stderr())
    if('fasta_file_path' %in% colnames(df)){
	df_complete = filter(df, !(is.na(fasta_file_path) | fasta_file_path == ''))
        df = filter(df, is.na(fasta_file_path) | fasta_file_path == '')
    } else {
        stop('Cannot find column: "fasta_file_path"')
    }
}

## Filtering based on user params & getting accessions
write('Filtering out genomes lacking NCBI genbank assembly accession', stderr())
col = unlist(args['column'])[1]
df = df[df[,col] != 'none',]
write(sprintf('Number of rows after filtering: %s', nrow(df)), stderr())
### just accessions
df_acc = df[,col]
df_acc = as.data.frame(df_acc)

# creating temp file of accessions
## Creating output directory
D = normalizePath(unlist(args['output'])[1])
dir.create(D, showWarnings = FALSE)
## writing table
F = file.path(D, 'accession.txt')
write(sprintf('Writing accessions to: %s', F), stderr())
write.table(df_acc, file=F, sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

# calling ncbi genome download
procs = as.character(unlist(args['procs'])[1])
retries = as.character(unlist(args['retries'])[1])
params = as.character(unlist(args['params'])[1])
database = as.character(unlist(args['database'])[1])
skip_bool = args['skip'][1]
if(skip_bool != TRUE){
    cmd = paste(c(exe, '-F', 'fasta', '-o', D, '-p', procs, '-r', retries,
                  '-A', F, '-s', database, params), collapse=' ')
    write(sprintf('Running cmd: %s', cmd), stderr())
    system(cmd)
} else {
    write('Skipping genome download', stderr())
}

# adding paths to genomes onto the table
## getting file paths
D2 = file.path(D, database)
fasta_files = list.files(D2, pattern='*.fna.gz', recursive=TRUE, full.names=TRUE)
n_files = as.character(length(fasta_files))
write(sprintf('Number of fasta files found: %s', n_files), stderr())
## Adding paths to input table
write('Adding file paths to the input table', stderr())
fasta_files = data.frame(accession = gsub('.+/', '', fasta_files),
                         fasta_file_path = fasta_files)
fasta_files$accession = gsub('(GCA_[0-9]+\\.[0-9]+)_.+', '\\1', fasta_files$accession)
df = left_join(df, fasta_files, by=setNames('accession', col))

# recombining tables (if --filter)
if(filter_bool == TRUE){
    df = rbind(df_complete, df)
}

# writing table
write.table(df, file=stdout(), sep='\t', row.names=FALSE, quote=FALSE)

# status
n_rows = as.character(nrow(df))
write(sprintf('Number of rows in the output: %s', n_rows), stderr())

n_missing = as.character(nrow(df[is.na(df$fasta_file_path),]))
write(sprintf('Number of rows with missing file paths: %s ', n_missing), stderr())



