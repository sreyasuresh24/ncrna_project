import pandas

os.makedirs('slurm',exist_ok=True)

sample_sheet = pandas.read_csv(config["sample_sheet"], sep="\t")

rule all:
    input:
        expand("fastqc/{library_name}_{read}_fastqc.html", library_name=sample_sheet["sra_id"], read=[1, 2]),
        expand("trim_galore/{library_name}_{read}_val_{read}.fq.gz", library_name=sample_sheet["sra_id"], read=[1, 2]),

rule fastqc:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sra_id == wc.sra_id]["fastq_1"],
    output:
        html="fastqc/{library_name}_1_fastqc.html"
    shell:
        r"""
        fastqc --outdir fastqc {input.fastq1}
        """

rule fastqc2:
    input:
        fastq2=lambda wc: sample_sheet[sample_sheet.sra_id == wc.sra_id]["fastq_2"],
    output:
        html="fastqc/{library_name}_2_fastqc.html"
    shell:
        r"""
        fastqc --outdir fastqc {input.fastq2}
        """
rule trim_galore:
    input:
        fastq1=lambda wc: sample_sheet[sample_sheet.sra_id == wc.sra_id]["fastq_1"], 
        fastq2=lambda wc: sample_sheet[sample_sheet.sra_id == wc.sra_id]["fastq_2"], 
    output:
        fastq1="trim_galore/{library_name}_1_val_1.fq.gz",  
        fastq2="trim_galore/{library_name}_2_val_2.fq.gz",
    shell:
        r"""
        trim_galore --phred64 --paired --output_dir trim_galore {input.fastq1} {input.fastq2}   
        """