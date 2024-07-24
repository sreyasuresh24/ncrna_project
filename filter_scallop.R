#!/usr/bin/env Rscript
library(data.table)
library(rtracklayer)
gtffile <- snakemake@input[['gtf']]
output <- snakemake@output[['gtf']]
sample_name <- snakemake@wildcards[['sample_name']]
gtf <- import(gtffile)
# Shorten transcript_id by truncation to a maximum of 50 characters
sn <- sample_name
err1 <- sub('.*ERR', 'ERR', sn)
err <- sub('\\..*|_.*', '', err1)
gtf$transcript_id <- paste0(err, "_", gtf$transcript_id)
gtf$gene_id <- paste0(err, "_", gtf$gene_id)
stopifnot(nchar(gtf$transcript_id) < 50)
dt <- as.data.table(gtf)

tx <- dt[type == 'exon', list(length=sum(width), n_exons=.N), transcript_id]
rpkm <- dt[type == 'transcript', list(transcript_id, RPKM)]
tx <- merge(tx, rpkm, by='transcript_id')
keep <- tx[length > 200 & RPKM > 5 & n_exons > 1]
outgtf <- export(gtf[gtf$transcript_id %in% keep$transcript_id], output)


