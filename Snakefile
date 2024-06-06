import pandas as pd
import os
import re
import sys

# Create the slurm directory if it doesn't exist
os.makedirs('slurm', exist_ok=True)

# Read the sample sheet
sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t")
sample_sheet['sample_name'] = [re.sub("_1.fastq.gz",'', os.path.basename(x)) for x in sample_sheet.fastq_1]
print(sample_sheet)

rule all:
    input:
        expand("fastqc/{sample_name}_{read}_fastqc.html", sample_name=sample_sheet["sample_name"], read=[1, 2]),
        expand("trim_galore/{sample_name}_{read}_val_{read}.fq.gz", sample_name=sample_sheet["sample_name"], read=[1, 2]),
        expand("hisat2_alignment/{sample_name}.bam", sample_name=sample_sheet["sample_name"]),

rule fastqc:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_1"],
    output:
        html="fastqc/{sample_name}_1_fastqc.html"
    shell:
        r"""
        fastqc --outdir fastqc {input.fastq1}
        """

rule fastqc2:
    input:
        fastq2=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_2"],
    output:
        html="fastqc/{sample_name}_2_fastqc.html"
    shell:
        r"""
        fastqc --outdir fastqc {input.fastq2}
        """

rule trim_galore:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_1"], 
        fastq2=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_2"], 
    output:
        fastq1="trim_galore/{sample_name}_1_val_1.fq.gz",  
        fastq2="trim_galore/{sample_name}_2_val_2.fq.gz",
    shell:
        r"""
        trim_galore --phred64 --paired --output_dir trim_galore {input.fastq1} {input.fastq2}   
        """

rule reference_genome:
    output:
        reference="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa",
        index="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index.1.ht2",
    params:
        index="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index",
    shell:
        r"""
        curl -L https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa.gz \
        | gunzip > {output.reference}
        hisat2-build {output.reference} {params.index}
        """
rule hisat_alignment:
    input:
        fastq1="trim_galore/{sample_name}_1_val_1.fq.gz",
        fastq2="trim_galore/{sample_name}_2_val_2.fq.gz",
        index="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index.1.ht2"
    output:
        bam="hisat2_alignment/{sample_name}.bam"
    threads: 4
    shell:
        r"""
        hisat2 -x reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index \
               -1 {input.fastq1} -2 {input.fastq2} \
               -p {threads} | samtools sort > {output.bam}
                """




        






    