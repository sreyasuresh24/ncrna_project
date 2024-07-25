library(data.table)
library(rtracklayer)

# File paths
bed_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.bed"
merged_gtf_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/merged_stringtie_scallop/merged.scallop.stringtie.5k.gtf"
output_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/filter/final_filtered_transcripts.gtf"

# Read the BED file
data <- fread(bed_file, skip = 1, select = c(1, 2, 3, 4))
setnames(data, names(data), c("tx_name", "Start", "End", "Attributes"))

# Extract score attribute
data[, score := sub(',.*', '', sub('.*score=', '', Attributes))]
data[, score := as.numeric(score)]

# List of coding transcripts (based on score > -Inf, meaning no cutoff)
exclude <- unique(data[score > -Inf, list(tx_name)])
# Read GFF file
gff_data <- fread(merged_gtf_file, sep="\t", header=FALSE, fill=TRUE, quote="")
setnames(gff_data, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"), 
         c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))

# Extract transcript IDs from the GFF data
gff_data[, tx_name := sub("ID=", "", sub(";.*", "", attributes))]

# Filter GFF data to exclude coding transcripts
filtered_gff <- gff_data[!tx_name %in% exclude$tx_name]

# Ensure GTF format conventions
filtered_gff[, attributes := sub("transcript_id", "gene_id", attributes)]
filtered_gff[, attributes := paste(attributes, paste0("transcript_id \"", tx_name, "\";"))]

# Write the filtered data to a GTF file
write.table(filtered_gff[, .(seqname, source, feature, start, end, score, strand, frame, attributes)], 
            file=output_file, 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
