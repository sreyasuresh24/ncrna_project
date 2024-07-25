library(data.table)

# File paths (make sure these paths are correctly defined)
bed_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.bed"
gff_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.gff3"
output_file <- "/mnt/data/project0014/Sreya/ncrna/output_results/filter/filter_transcript.gtf"

# Read BED file
data <- fread(bed_file, skip = 1, select = c(1, 2, 3, 4))
setnames(data, names(data), c("tx_name", "Start", "End", "Attributes"))

# Extract score from Attributes and convert to numeric
data[, score := as.numeric(sub(',.*', '', sub('.*score=', '', Attributes)))]

# Identify transcripts to exclude based on score
exclude <- unique(data[score > -Inf, list(tx_name)])

# Read the GFF file
gff_data <- fread(gff_file, header = FALSE, sep = "\t", quote = "", fill = TRUE)

# Set column names
setnames(gff_data, c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

# Extract transcript IDs from the GFF file (assuming IDs are in the attributes column)
gff_data[, tx_name := sub(".*ID=([^;]+);.*", "\\1", attributes)]

# Filter out transcripts present in the exclude list
filtered_gff <- gff_data[!tx_name %in% exclude$tx_name]

# Check if output directory exists and create it if it doesn't
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the filtered GFF file in GTF format (ensure to use GTF format specifics if needed)
fwrite(filtered_gff, output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
 
