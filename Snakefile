import pandas as pd
import os
import re
import sys

os.makedirs('slurm', exist_ok=True)
sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t", comment='@')
sample_sheet['sample_name'] = [re.sub("_1.fastq.gz",'', os.path.basename(x)) for x in sample_sheet.fastq_1]



rule all:
    input:
        expand("fastqc/{sample_name}_{read}_fastqc.html", sample_name=sample_sheet["sample_name"], read=[1, 2]),
        expand("trim_galore/{sample_name}_{read}_val_{read}.fq.gz", sample_name=sample_sheet["sample_name"], read=[1, 2]),
        expand("hisat2_alignment/{sample_name}.bam", sample_name=sample_sheet["sample_name"]),
        expand("stringtie_assembly/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        expand("scallop_assembly/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        expand("filter_stringtie/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        expand("filter_scallop/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        "merged_stringtie_scallop/merged.scallop.stringtie.5k.gtf",
        "transcripts/transcripts.5k.withgffread.fa",
        "transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder_dir/longest_orfs.pep",
        "6157.faa",
        "6157.faa.dmnd",
        "blast.outfmt6",
        "Pfam-A.hmm",
        "pfam.domtblout",
        "transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.gff3",
        "transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.pep",
        "transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.cds",
        "repeatmasker/transcripts.5k.withgffread.fa.cat.gz",
        "repeatmasker/transcripts.5k.withgffread.fa.out",
        "repeatmasker/transcripts.5k.withgffread.fa.tbl",
        "repeatmasker/transcripts.5k.withgffread.fa.masked",
        "repeatmasker/transcripts.5k.withgffread.fa.out.gff",
        "filter/filter_transcript.gtf"

rule fastqc:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_1"]  
    output:
        html="fastqc/{sample_name}_1_fastqc.html"
    shell:
        "fastqc --outdir fastqc {input.fastq1}"

rule fastqc2:
    input:
        fastq2=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_2"]
    output:
        html="fastqc/{sample_name}_2_fastqc.html"
    shell:
        "fastqc --outdir fastqc {input.fastq2}"

rule trim_galore:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_1"],
        fastq2=lambda wc: sample_sheet[sample_sheet.sample_name == wc.sample_name]["fastq_2"]
    output:
        fastq1="trim_galore/{sample_name}_1_val_1.fq.gz",
        fastq2="trim_galore/{sample_name}_2_val_2.fq.gz"
    shell:
        "trim_galore --phred33 --paired --output_dir trim_galore {input.fastq1} {input.fastq2}"

rule reference_genome:
    output:
        reference="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa",
        index="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index.1.ht2"
    params:
        index="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index"
    shell:
        """
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
    resources:
        mem='20G'
    threads: 4
    shell:
        """
        hisat2 -x reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked_index \
               -1 {input.fastq1} -2 {input.fastq2} \
               -p {threads} | samtools sort > {output.bam}
        """
rule reference_annotations:
    output:
        reference_annotation='reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3'
    shell:
        r"""
        curl -L https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz \
        """
    
rule stringtie_assembly:
    input:
        bam="hisat2_alignment/{sample_name}.bam",
        gff="reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
    output:
        gtf="stringtie_assembly/{sample_name}.gtf"
    resources:
        cpus=16,
        mem='20G'
    shell:
        r"""
        stringtie {input.bam} -o {output.gtf} -p {resources.cpus} -G {input.gff}
        """

rule scallop_assembly:
    input:
        bam="hisat2_alignment/{sample_name}.bam",
        gff="reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
    output:
        gtf="scallop_assembly/{sample_name}.gtf"
    resources:
        cpus=16,
        mem='40G'
    shell:
        r"""
        scallop -i {input.bam} -o {output.gtf} --gtf {input.gff} --threads {resources.cpus}
        """

rule filter_stringtie:
    input:
        gtf="stringtie_assembly/{sample_name}.gtf"
    output:
        gtf="filter_stringtie/{sample_name}.gtf"
    script:
        "../scripts/filter_stringtie.R"

rule filter_scallop:
    input:
        gtf="scallop_assembly/{sample_name}.gtf"
    output:
        gtf="filter_scallop/{sample_name}.gtf"
    script:
        "../scripts/filter_scallop.R"

rule merge_stringtie_scallop:
    input:
        stringtie_gtf = expand("filter_stringtie/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        scallop_gtf = expand("filter_scallop/{sample_name}.gtf", sample_name=sample_sheet["sample_name"]),
        reference = "reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa",
    output:
        merged_gtf="merged_stringtie_scallop/merged.scallop.stringtie.5k.gtf",
    shell:
        r"""
        gffread {input.stringtie_gtf} {input.scallop_gtf} -N -M -F -K -Q -Z -i 5000 -T -g {input.reference} -o {output.merged_gtf}
        """

rule extract_fasta:
    input:
        merged_gtf="merged_stringtie_scallop/merged.scallop.stringtie.5k.gtf",
        genome="reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa",
    output:
        fasta="transcripts/transcripts.5k.withgffread.fa",
    shell:
        r"""
        gffread -w {output.fasta} -g {input.genome} {input.merged_gtf}
        """

rule transdecoder_longorfs:
    input:
        "transcripts/transcripts.5k.withgffread.fa"
    output:
        "transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder_dir/longest_orfs.pep"
    shell:
        """
        mkdir -p transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir
        TransDecoder.LongOrfs -m 30 -t {input} --output_dir transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir
        """
rule download_genome_proteins:
    output:
        "6157.faa"
    shell:
        r"""
        datasets download genome taxon 6157 --dehydrated --include protein --filename 6157-genome-protein.zip
        unzip 6157-genome-protein.zip -d 6157-genome-protein
        datasets rehydrate --directory 6157-genome-protein
        cat 6157-genome-protein/ncbi_dataset/data/*/protein.faa > 6157.faa
        """

rule diamond_makedb:
    input:
        "6157.faa"
    output:
        "6157.faa.dmnd"
    threads: 8  
    shell:
        r"""
        diamond makedb --in {input} -d {output}
        """
rule diamond_blastp:
    input:
        query="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder_dir/longest_orfs.pep",
        database="6157.faa.dmnd"
    output:
        "blast.outfmt6"
    threads: 8  
    shell:
        r"""
        diamond blastp --query-cover 90 --id 80 --query {input.query} --outfmt 6 --max-target-seqs 1 --evalue 1e-5 --threads {threads} -o {output} --db {input.database}
        """

rule download_pfam_database:
    output:
        "Pfam-A.hmm"
    shell:
        r"""
        curl -L http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz | gunzip > {output}
        """

rule hmmer_search:
    input:
        pfam_db="Pfam-A.hmm",
        query="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder_dir/longest_orfs.pep"
    output:
        "pfam.domtblout"
    threads: 8  
    shell:
        r"""
        hmmsearch --cpu {threads} -E 1e-5 --domtblout {output} {input.pfam_db} {input.query} > /dev/null
        """

rule transdecoder_predict:
    input:
        transcripts="transcripts/transcripts.5k.withgffread.fa",
        blast_hits="blast.outfmt6",
        pfam_hits="pfam.domtblout"
    output:
        gff3="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.gff3",
        pep="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.pep",
        cds="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.cds"
    params:
        outdir="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir"
    resources:
        mem='10G'
    shell:
        """
        TransDecoder.Predict -t {input.transcripts} --output_dir {params.outdir} --retain_blastp_hits {input.blast_hits} --retain_pfam_hits {input.pfam_hits}
        """

rule repeatmasker:
    input:
        "transcripts/transcripts.5k.withgffread.fa"
    output:
        cat_gz="repeatmasker/transcripts.5k.withgffread.fa.cat.gz",
        out="repeatmasker/transcripts.5k.withgffread.fa.out",
        tbl="repeatmasker/transcripts.5k.withgffread.fa.tbl",
        masked="repeatmasker/transcripts.5k.withgffread.fa.masked",
        out_gff="repeatmasker/transcripts.5k.withgffread.fa.out.gff"
    shell:
        """
        mkdir -p repeatmasker
        RepeatMasker -pa 24 -gff -dir repeatmasker -species Platyhelminthes {input}
        """

rule filter_gff:
    input:
        bed_file="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.bed",
        gff_file="transdecoder/transcripts.5k.withgffread.fa.transdecoder_dir/transcripts.5k.withgffread.fa.transdecoder.gff3"
    output:
        filtered_gff="filter/filter_transcript.gtf"
    script:
        "filter.R"



    