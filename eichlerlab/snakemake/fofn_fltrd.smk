# This pipeline filters out low QV reads in the cell and outputs a subsetted fastq
import pandas as pd
import numpy as np

manifest_fp = 'fofn_fltrd.tsv'
manifest_df = pd.read_csv(manifest_fp, sep='\t', index_col='cell_name')

wildcard_constraints:
    sample='|'.join(manifest_df['sample']),
    cell_name='|'.join(manifest_df.index),


def find_fastq(which_fp='fastq'):
    def inner(wildcards):
        fofn_fp = manifest_df.at[wildcards.cell_name, 'fofn']
        fofn_df = pd.read_csv(fofn_fp, names=['fastq'])

        fastq_fp = fofn_df.loc[fofn_df['fastq'].str.contains(wildcards.cell_name), 'fastq'].tolist()[0]
        if which_fp == 'fastq':
            return fastq_fp
        else:
            return fastq_fp + '.fai'
    return inner

def get_fofn(wildcards):
    return manifest_df.at[wildcards.cell_name, 'fofn']

def get_qv(wildcards):
    return manifest_df.at[wildcards.cell_name, 'qv_fp']

rule all:
    input:
        expand('{sample}_{cell_name}.fastq.fofn', zip, sample=manifest_df['sample'], cell_name=manifest_df.index)

rule extract_lowQ_reads:
    input:
        yak_out_txt = get_qv,
        fastq_fai = find_fastq(which_fp='fai')
    output:
        lowQ_reads = '{cell_name}_lowQV-reads.txt',
        target_reads = '{cell_name}_target-reads.txt',
        new_fai = '{cell_name}-subset.fastq.gz.fai'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    run:
        yak_out_df = pd.read_csv(input.yak_out_txt,sep='\t',comment='C',header=None,usecols=[1, 5],names=['read_name', 'qv']).dropna()
        low_quality = np.mean(yak_out_df['qv']) - np.std(yak_out_df['qv'])

        yak_out_df = yak_out_df[yak_out_df['qv'] < low_quality]
        yak_out_df.to_csv(output.lowQ_reads, index=False, header=False, sep='\t')

        fastq_fai_df = pd.read_csv(input.fastq_fai, sep='\t', header=None, names=['read_name', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'])
        fastq_fai_df = fastq_fai_df[~fastq_fai_df['read_name'].isin(yak_out_df['read_name'])]
        fastq_fai_df.to_csv(output.target_reads,header=False,index=False,columns=['read_name'], sep='\t')
        fastq_fai_df.to_csv(output.new_fai,header=False,index=False, sep='\t')

parts = 15
scattergather:
    sg_parts = parts

rule split_fastq:
    input:
        target_reads=rules.extract_lowQ_reads.output.target_reads,
    output:
        splitted = scatter.sg_parts('temp/{{cell_name}}_target-reads_{scatteritem}.txt')
    threads: 1
    resources:
        mem = 1,
        hrs=12
    run:
        df = pd.read_csv(input.target_reads, header=None)

        dfs = np.array_split(df, parts)
        for idx, entry in enumerate(dfs):
            entry.to_csv(output.splitted[idx], index=False, header=False)


rule subset_fastq:
    input:
        target_reads=rules.split_fastq.output.splitted,
        fastq=find_fastq(which_fp='fastq')
    output:
        subsetted_fastq = 'temp/{cell_name}-subset_{scatteritem}.fastq.gz'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    container: 'docker://quay.io/biocontainers/seqtk:1.3--h7132678_4'
    shell:
        '''
        seqtk subseq {input.fastq} {input.target_reads} | gzip -c - > {output.subsetted_fastq}
        '''

rule gather_fastq:
    input:
        subset_fastqs = gather.sg_parts('temp/{{cell_name}}-subset_{scatteritem}.fastq.gz')
    output:
        new_fastq = '{cell_name}-subset.fastq.gz'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        zcat {input.subset_fastqs} | gzip -c - > {output.new_fastq}
        '''

rule make_new_fofn:
    input:
        fofn = get_fofn,
        fastq = rules.gather_fastq.output.new_fastq
    output:
        new_fofn = '{sample}_{cell_name}.fastq.fofn'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        cat <(grep -v "{wildcards.cell_name}" {input.fofn}) <( readlink -f {input.fastq} ) > {output.new_fofn}
        '''

onsuccess: shell('rm -rf temp/')