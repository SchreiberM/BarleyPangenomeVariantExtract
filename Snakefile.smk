#!/bin/bash -e

import pandas as pd
import os
 
#dir = "workflow/logs/"
if not os.path.exists("workflow/logs/"):
    os.mkdir("workflow/logs/")

configfile: "config.yaml"

genomeseq=config["GenomeSequence"]
referenceAlignment=config["ReferenceToAlignTo"]

genomes = pd.read_table("workflow/blastLocations/PangenomeBlast_Sub.tsv").set_index("genotype", drop=False)
sequences = pd.read_table("workflow/blastLocations/PangenomeSequences_Sub.tsv").set_index("genotype", drop=False)


def blastDatabase(wildcards):
    return expand(
        "{reference}", reference=genomes["path"][wildcards.sample])

def genomeRef(wildcards):
    return expand(
        "{reference}", reference=sequences["path"][wildcards.sample])

def get_genomeSequenceInput(wildcards):
    return config["GenomeSequence"][wildcards.genomeseq]


rule all:
    input:
        expand( "workflow/{genomeseq}/report/AllVariantCalls.{referenceAlignment}.vcf", genomeseq=config["GenomeSequence"],referenceAlignment=config["ReferenceToAlignTo"] ),
        expand( "workflow/{genomeseq}/report/Figure1.pdf",  genomeseq=config["GenomeSequence"])


rule blast_nucleotide:
    input:
        query = get_genomeSequenceInput
    output:
        "workflow/{genomeseq}/blastResult/{sample}.blast.txt"
    log:
        "workflow/logs/{genomeseq}/blastsearch/{sample}.blast.log"
    threads:
        1
    conda:
        "workflow/envs/panAnalysis.yml"
    params:
        format="'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand'",
        percent=config["PercentIdentity"],
        coverage=config["PercentCoverage"],
        blastdb=blastDatabase
    shell:
        "blastn -query {input.query} -outfmt {params.format} -perc_identity {params.percent} -qcov_hsp_perc {params.coverage} -db {params.blastdb} -num_threads {threads} -out {output}"


rule extract_coordinates:
    input: 
        files= expand( "workflow/{{genomeseq}}/blastResult/{sample}.blast.txt", sample=genomes["genotype"] ),
        query = get_genomeSequenceInput
    output:
        Bedfiles=directory("workflow/{genomeseq}/bedFiles/"),
        figure=report( "workflow/{genomeseq}/report/Figure1.pdf" )
    log:
        "workflow/logs/{genomeseq}/ExtractCoordinates.log"
    threads:
        1
    conda:
        "workflow/envs/panAnalysis.yml"
    params:
        distance=config["DistanceUpDownstream"],
        filterIdentity=config["PercentIdentityFilter"],
        filterHsp=config["PercentCoverageFilter"],
        coverage=config["PercentCoverage"],
    priority:
        50
    shell:
        """
        Rscript workflow/rscripts/plotCheckpoint.r {output.figure} {params.distance} {output.Bedfiles} {input.query} {params.filterIdentity} {params.filterHsp} {params.coverage}
        """

rule extract_Sequences:
    input: 
        genomeFile=genomeRef,
        figure="workflow/{genomeseq}/report/Figure1.pdf"
    output:
        finalSequences="workflow/{genomeseq}/sequences/{sample}.fasta"
    log:
        "workflow/logs/{genomeseq}/{sample}.bedfile.log"
    threads:
        1
    conda:
        "workflow/envs/panAnalysis.yml"
    params:
        "workflow/{genomeseq}/bedFiles/{sample}.bed",
    shell:
        """
        bedtools getfasta -nameOnly -s -fo {output.finalSequences} -fi {input.genomeFile} -bed {params}
        """

rule minimap2_index:
    input:
        target="workflow/{genomeseq}/sequences/{referenceAlignment}.fasta"
    output:
        "workflow/{genomeseq}/minimap/{referenceAlignment}.mmi"
    log:
        "workflow/logs/{genomeseq}/minimap2_index/{referenceAlignment}.log"
    params:
        extra=""  # optional additional args
    threads: 4
    wrapper:
        "v8.1.1/bio/minimap2/index"


rule minimap2_bam_sorted:
    input:
        target="workflow/{genomeseq}/minimap/{referenceAlignment}.mmi",  # can be either genome index or genome fasta
        query="workflow/{genomeseq}/sequences/{sample}.fasta"
    output:
        "workflow/{genomeseq}/aligned/{sample}.{referenceAlignment}.aln.sorted.bam"
    log:
        "workflow/logs/{genomeseq}/minimap2/{sample}.{referenceAlignment}.log"
    params:
        extra=r"""-R '@RG\tID:{sample}\tSM:{sample}' -x asm20 -N 0""",  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 2
    wrapper:
        "v8.1.1/bio/minimap2/aligner"

rule mergeBam:
    input:
        expand("workflow/{{genomeseq}}/aligned/{sample}.{referenceAlignment}.aln.sorted.bam", sample=genomes["genotype"], referenceAlignment=config["ReferenceToAlignTo"])
    output:
        "workflow/{genomeseq}/variant/merged.bam"
    log:
        "workflow/logs/{genomeseq}/mergeBam.log"
    threads: 2
    conda:
        "workflow/envs/panAnalysis.yml"
    shell:
        """
        samtools merge {output} {input}
        samtools index {output}
        """

rule nvc_variants:
    input:
        bamFile="workflow/{genomeseq}/variant/merged.bam",
        target="workflow/{genomeseq}/sequences/{referenceAlignment}.fasta"
    output:
        "workflow/{genomeseq}/report/AllVariantCalls.{referenceAlignment}.vcf"
    log:
        "workflow/logs/{genomeseq}/{referenceAlignment}.nvc.log"
    threads: 4
    conda:
        "workflow/envs/panAnalysis.yml"
    shell:
        """
        samtools faidx {input.target}
        naive_variant_caller.py -b {input.bamFile} -o {output} -r {input.target} --variants_only --min_base_quality=0 --min_mapping_quality=0 --ploidy=1
        """
