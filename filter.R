library(data.table)

# Input and output file paths provided by Snakemake
bed_file <- snakemake@input[["bed_file"]]
gff_file <- snakemake@input[["gff_file"]]
output_file <- snakemake@output[["filtered_gff"]]

# Read BED file
bed_data <- fread(bed_file, skip=1, select=c(1,2,3,4))
setnames(bed_data, c("V1", "V2", "V3", "V4"), c("tx_name", "Start", "End", "Attributes"))
bed_data[, score := as.numeric(sub(',.*', '', sub('.*score=', '', Attributes)))]

# List of coding transcripts to exclude
exclude <- unique(bed_data[score > -Inf, .(tx_name)])

# Read GFF file
gff_data <- fread(gff_file, sep="\t", header=FALSE, fill=TRUE, quote="")
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

