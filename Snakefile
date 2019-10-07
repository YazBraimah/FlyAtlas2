"""
Author: Y. Ahmed-Braimah
--- RNA-seq analysis workflow for FlyAtlas2 mRNAs and miRNAs.

"""

import json
import os
import re
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']
GTF = config['GTF']
TRX = config['TRX']
# MIR = config['MIR']
INDEX = config['INDEX']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
WORK_DIR = config['WORK_DIR']
HOME_DIR = config['HOME_DIR']  # the "launch_snakemake.sh" and "config.yml" files should be here

## set the usr and job environments for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')

# Samples and their corresponding filenames.
# paired-end mRNA
peFILESmRNA = json.load(open(config['PE_SAMPLES_mRNA'])) 
peSAMPLESmRNA = sorted(peFILESmRNA.keys())                  
# single-end microRNA
# seFILESmicroRNA = json.load(open(config['SE_SAMPLES_microRNA'])) 
# seSAMPLESmicroRNA = sorted(seFILESmicroRNA.keys())
# # paired-end microRNA
# peFILESmicroRNA = json.load(open(config['PE_SAMPLES_microRNA'])) 
# peSAMPLESmicroRNA = sorted(peFILESmicroRNA.keys())     

# combined single and paired microRNA samples
# FILES = json.load(open(config['SAMPLES_JSON']))
# combinedSamMrna = [seSAMPLESmicroRNA, peSAMPLESmicroRNA]
# SAMPLESmicroRNA = [y for x in combinedSamMrna for y in x]

# combine ALL samples
# combinedSam = [peSAMPLESmRNA, seSAMPLESmicroRNA, peSAMPLESmicroRNA]
# SAMPLES = [y for x in combinedSam for y in x]

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

##--------------------------------------------------------------------------------------##
#
# _____ _             _               _               _       
#|  ___(_)_ __   __ _| |   ___  _   _| |_ _ __  _   _| |_ ___ 
#| |_  | | '_ \ / _` | |  / _ \| | | | __| '_ \| | | | __/ __|
#|  _| | | | | | (_| | | | (_) | |_| | |_| |_) | |_| | |_\__ \
#|_|   |_|_| |_|\__,_|_|  \___/ \__,_|\__| .__/ \__,_|\__|___/
#                                        |_|                  
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all: 
    input:
        join(OUT_DIR, 'ballgown', 'gene_counts.csv'), 
        join(OUT_DIR, 'ballgown', 'transcript_counts.csv'),
        # join(OUT_DIR, 'eXpress_microRNA', 'express.gene.TMM.EXPR.matrix'),
        # join(OUT_DIR, 'eXpress_mRNA', 'express.gene.TMM.EXPR.matrix'),
        # join(OUT_DIR, 'MultiQC_microRNA', 'multiqc_report.html'),
        join(OUT_DIR, 'MultiQC_mRNA', 'multiqc_report.html'),
        join(OUT_DIR, 'StringTie', 'gffcmp.annotated.gtf')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
#  ___   ____ 
# / _ \ / ___|
#| | | | |    
#| |_| | |___ 
# \__\_\\____|
#             
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# ## Rule to check raw microRNA SE read quality
# rule fastqcSE_microRNA:
#     input:
#         r1 = lambda wildcards: seFILESmicroRNA[wildcards.sample]['R1']
#     output:
#         r1 = join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.R1_fastqc.html')
#     log:
#         join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.fastQC_se.log')
#     benchmark:
#         join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.fastQC_se.benchmark.tsv')
#     message: 
#         """--- Checking read quality of SE sample "{wildcards.sample}" with FastQC """
#     run:
#         if not os.path.exists(join(OUT_DIR, 'fastQC_mmicroRNA')):
#             os.makedirs(join(OUT_DIR, 'fastQC_microRNA'))

#         shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
#                 ' && cp {input.r1} ' + join(WORK_DIR, USER, JOB_ID) +
#                 ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
#                 ' && fastqc {wildcards.sample}.R1.fq.gz' 
#                 ' > {log} 2>&1'
#                 ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC_microRNA'))
#         shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# ## Rule to check raw microRNA PE read quality
# rule fastqcPE_microRNA:
#     input:
#         r1 = lambda wildcards: peFILESmicroRNA[wildcards.sample]['R1'],
#         r2 = lambda wildcards: peFILESmicroRNA[wildcards.sample]['R2']
#     output:
#         r1 = join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.R1_fastqc.html'),
#         r2 = join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.R2_fastqc.html')
#     log:
#         join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.fastQC_init_pe.log')
#     benchmark:
#         join(OUT_DIR, 'fastQC_microRNA', '{sample}' + '.fastQC_init_pe.benchmark.tsv')
#     message: 
#         """--- Checking read quality of PE sample "{wildcards.sample}" with FastQC """
#     run:
#         if not os.path.exists(join(OUT_DIR, 'fastQC_microRNA')):
#             os.makedirs(join(OUT_DIR, 'fastQC_microRNA'))

#         shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
#                 ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
#                 ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
#                 ' && fastqc {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz' 
#                 ' > {log} 2>&1'
#                 ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC_microRNA'))
#         shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule fastqcPE_mRNA:
    input:
        r1 = lambda wildcards: peFILESmRNA[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILESmRNA[wildcards.sample]['R2']
    output:
        r1 = join(OUT_DIR, 'fastQC_mRNA', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC_mRNA', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC_mRNA', '{sample}' + '.fastQC_init_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC_mRNA', '{sample}' + '.fastQC_init_pe.benchmark.tsv')
    message: 
        """--- Checking read quality of PE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC_mRNA')):
            os.makedirs(join(OUT_DIR, 'fastQC_mRNA'))

        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC_mRNA'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
#
# _____                   _      ____  
#|_   _|   ___  _____  __| | ___|___ \ 
#  | || | | \ \/ / _ \/ _` |/ _ \ __) |
#  | || |_| |>  <  __/ (_| | (_) / __/ 
#  |_| \__,_/_/\_\___|\__,_|\___/_____|
#                                      
#
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to map PE reads with HISAT2
rule hisat2:
    input:
        r1 = lambda wildcards: peFILESmRNA[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILESmRNA[wildcards.sample]['R2']
    output:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}.csorted.hisat2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hisat2_{sample}.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_pe.benchmark.tsv')
    message: 
        """--- Mapping PE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(INDEX, '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID))
        shell('cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (hisat2'
                ' -p 8'
                ' --no-unal'
                ' --trim5 10'
                ' --dta'
                ' -x ' + join(rstrip(os.path.basename(DNA), '.fasta') + '_tran') +
                ' -1 {wildcards.sample}.R1.fq.gz'
                ' -2 {wildcards.sample}.R2.fq.gz)'
                ' 2>{log}'
                ' | samtools sort -@ 8 -o csorted.hisat2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.hisat2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.hisat2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to assemble transcripts with StringTie
rule stringtie:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.hisat2.bam')
    output:
        asmbly = join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.log')
    benchmark:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.benchmark.tsv')
    message: 
        """--- Assembling transcripts for sample {wildcards.sample} with StringTie """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' {wildcards.sample}.csorted.hisat2.bam'
                ' -p 8'
                ' -G ' + os.path.basename(GTF) +
                ' -o {wildcards.sample}.gtf'
                ' -l {wildcards.sample} > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.gtf') + ' ' + join(OUT_DIR, 'StringTie', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule merge_assemblies:
    input:
        assemblies = expand(join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf'), sample = peSAMPLESmRNA)
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'stringtie_merged.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', 'st_mrg.log')
    benchmark:
        join(OUT_DIR, 'StringTie', 'st_mrg.benchmark.tsv')
    message: 
        """--- Merging StringTie transcripts """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && ls -1 ' + join(OUT_DIR) + '/StringTie/*/*.gtf > ' + join(OUT_DIR, 'StringTie', 'assemblies.txt') +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' --merge'
                ' -p 4'
                ' -G ' + os.path.basename(GTF) +
                ' -o stringtie_merged.gtf ' + join(OUT_DIR, 'StringTie', 'assemblies.txt'))
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'stringtie_merged.gtf') + ' ' + join(OUT_DIR, 'StringTie'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule compare_gtf:
    input:
        STasmbly = rules.merge_assemblies.output.asmbly
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'gffcmp.annotated.gtf')
    params:
        gtf = GTF
    message: 
        """--- Comparing StringTie merged assembly with reference GTF """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} {input.STasmbly} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && gffcompare'
                ' -G'
                ' -r ' + os.path.basename(GTF) +
                ' stringtie_merged.gtf')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'gffcmp.*') + ' ' + join(OUT_DIR, 'StringTie'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to measure transcript abundances with Stringtie
rule abundances:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.hisat2.bam'),
        mrgd = rules.merge_assemblies.output.asmbly
    output:
        abundance = join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf')
    log:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.log')
    benchmark:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.benchmark.tsv')
    message: 
        """--- Estimating transcript abundances for sample {wildcards.sample} with StringTie"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {input.mrgd} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' -e -B -p 8'
                ' -G stringtie_merged.gtf' 
                ' -o {wildcards.sample}_abundance.gtf'
                ' {wildcards.sample}.csorted.hisat2.bam'
                ' > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}_abundance.gtf') + ' ' + join(OUT_DIR, 'ballgown', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# Rule to combine abundance counts for downstream analysis
rule collate_counts:
    input:
        abundances = expand(join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf'), sample = peSAMPLESmRNA)
    output:
        geneCounts = join(OUT_DIR, 'ballgown', 'gene_counts.csv'),
        transcriptCounts = join(OUT_DIR, 'ballgown', 'transcript_counts.csv')
    message: 
        """--- Outputting count matrices """
    run:
        shell('prepDE.py'
                ' -i ' + join(OUT_DIR, 'ballgown') + 
                ' -g ' + join(OUT_DIR, 'ballgown', 'gene_counts.csv') +
                ' -t ' + join(OUT_DIR, 'ballgown', 'transcript_counts.csv'))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
#     __  __                        
#  ___\ \/ /_ __  _ __ ___  ___ ___ 
# / _ \\  /| '_ \| '__/ _ \/ __/ __|
#|  __//  \| |_) | | |  __/\__ \__ \
# \___/_/\_\ .__/|_|  \___||___/___/
#          |_|                      
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


# Rule for making a cDNA file from the new GTF file and the DNA, then prep index for bowtie2.
rule index_mRNA:
    input:
        trx = TRX
    output:
        geneTrans = join(OUT_DIR, 'eXpress_mRNA', 'transcriptome', join(rstrip(os.path.basename(TRX), '.fasta') + '.gene_trans_map'))
    log:
        trans_bt2 = join(OUT_DIR, 'eXpress_mRNA', 'transcriptome', 'logs', 'trans_bt2.log')
    benchmark:
        join(OUT_DIR, 'eXpress_mRNA', 'transcriptome', 'logs', 'gffread.benchmark.tsv')
    message: 
        "--- Building bowtie2 transcriptome index for gffread transcripts."
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cp {input.trx} ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cd ' + join(WORK_DIR, USER, JOB_ID) +
        # Extract the FASTA header from the cDNA file and make into
        # trans_map file.
              ' && grep ">" ' + os.path.basename(TRX) + ' | sed "s/>//g" | sed "s/ .*parent=/ /g" | sed "s/;.*//g" | awk \'{{print $2"\t"$1}}\' | sort -u > ' + join(rstrip(os.path.basename(TRX), '.fasta') + '.gene_trans_map') +
        # And finally make the index files.
              ' && align_and_estimate_abundance.pl' 
              ' --transcripts ' + os.path.basename(TRX) +
              ' --gene_trans_map ' + join(rstrip(os.path.basename(TRX), '.fasta') + '.gene_trans_map') +
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --prep_reference'
              ' --output_dir transcriptome/'
              ' > {log.trans_bt2} 2>&1')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, join(rstrip(os.path.basename(TRX), '.fasta') + '*')) + ' ' + join(OUT_DIR, 'eXpress_mRNA', 'transcriptome'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# Rule for mapping SE reads to the new transcriptome file with bowtie2 and quantifying abundance with eXpress
rule eXpress:
    input:
        r1 = lambda wildcards: peFILESmRNA[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILESmRNA[wildcards.sample]['R2'],
        geneTrans = rules.index_mRNA.output.geneTrans
    output:
        results = join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'results.xprs.genes'),
        aln_stats = join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'bowtie2_{sample}.log')
    log:
        join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'logs', 'bowtie2_eXpress.log')
    benchmark:
        join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'logs', 'eXpress.benchmark.tsv')
    message: 
        """--- Mapping "{wildcards.sample}" reads to transcriptome with bowtie2 and quantifying abundance with eXpress."""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cp ' + join(OUT_DIR, 'eXpress_mRNA', 'transcriptome', join(rstrip(os.path.basename(TRX), '.fasta') + '*')) + ' ' + join(WORK_DIR, USER, JOB_ID))
        shell('cd ' + join(WORK_DIR, USER, JOB_ID) +
              ' && align_and_estimate_abundance.pl' 
              ' --transcripts ' + os.path.basename(TRX) +
              ' --seqType fq'
              ' --left {wildcards.sample}.R1.fq.gz'
              ' --right {wildcards.sample}.R2.fq.gz'
              ' --gene_trans_map ' + join(rstrip(os.path.basename(TRX), '.fasta') + '.gene_trans_map') +
              ' --thread_count 4'  
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --output_dir {wildcards.sample}'
              ' --bowtie2_eXpress "--trim5 10"'
              ' > {log} 2>&1')
        shell('awk \'/CMD: set/{{f=1;next}} /CMD: touch/{{f=0}} f\' {log} | sed "/bam/d" > {output.aln_stats}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}') + '/* ' + join(OUT_DIR, 'eXpress_mRNA', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))   

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule merge_abundance:
    input:
        quants = expand(join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'results.xprs.genes'), sample = peSAMPLESmRNA),
        geneTrans = join(OUT_DIR, 'eXpress_mRNA', 'transcriptome', join(rstrip(os.path.basename(TRX), '.fasta') + '.gene_trans_map'))
    output:
        abundances = join(OUT_DIR, 'eXpress_mRNA', 'express.gene.TMM.EXPR.matrix'),
        samplesList = join(OUT_DIR, 'eXpress_mRNA', 'isoforms.samples.list')
    log:
        join(OUT_DIR, 'eXpress_mRNA', 'abnd_merge.log')
    benchmark:
        join(OUT_DIR, 'eXpress_mRNA', 'abnd_merge.benchmark.tsv')        
    message: 
        "--- Merging eXpress outputs from all samples"
    run:
        shell('ls -1 ' + join(OUT_DIR, 'eXpress_mRNA', '*', 'results.xprs') + ' > ' + join(OUT_DIR, 'eXpress_mRNA', 'isoforms.samples.list'))
        shell('cd ' + join(OUT_DIR, 'eXpress_mRNA') +
                ' && abundance_estimates_to_matrix.pl'
                ' --est_method eXpress'
                ' --gene_trans_map {input.geneTrans}'
                ' --name_sample_by_basedir'
                ' --out_prefix express'
                ' --quant_files ' + join(OUT_DIR, 'eXpress_mRNA', 'isoforms.samples.list') +
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
# __  __       _ _   _  ___   ____ 
#|  \/  |_   _| | |_(_)/ _ \ / ___|
#| |\/| | | | | | __| | | | | |    
#| |  | | |_| | | |_| | |_| | |___ 
#|_|  |_|\__,_|_|\__|_|\__\_\\____|
#                                  
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.hisat2.bam'), sample = peSAMPLESmRNA),
        # expand(join(OUT_DIR, 'eXpress_mRNA', '{sample}', 'results.xprs.genes'), sample = peSAMPLESmRNA),
        expand(join(OUT_DIR, 'fastQC_mRNA', '{sample}' + '.R1_fastqc.html'), sample = peSAMPLESmRNA)
        
    output:
        join(OUT_DIR, 'MultiQC_mRNA', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC_mRNA', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC_mRNA', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/HISAT-2/*/*log > ' + join(OUT_DIR, 'MultiQC_mRNA', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC_mRNA/*fastqc.zip >> ' + join(OUT_DIR, 'MultiQC_mRNA', 'summary_files.txt'))
        # shell('ls -1 ' + join(OUT_DIR) + '/eXpress_mRNA/*/bowtie2_*.log >> ' + join(OUT_DIR, 'MultiQC_mRNA', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC_mRNA') + ' -l ' + join(OUT_DIR, 'MultiQC_mRNA', 'summary_files.txt') +
                ' > {log} 2>&1')


