# sequence of tandem repeat with some flanking sequence
# periodicity (i.e. motif length)
# establish abundance of each available k-mer that has the same length as motif
# order kmers by decreasing order of abundance

import pandas as pd

configfile: 'config.yaml'

manifest = config.get("manifest", "manifest.tab") # sample\thap1\thap2
manifest_df = pd.read_csv(manifest, sep='\t', index_col='sample')

target_region = config.get('target_region', 'target.bed')
minimap_params = config.get('minimap_params', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B5')
motif = config['motif']
ref = config['ref']
ref_path = config['ref_path']

wildcard_constraints:
    sample = '|'.join(manifest_df.index),
    hap = 'hap1|hap2'

def get_fasta(wildcards):
    return manifest_df.at[wildcards.sample, wildcards.hap]

rule all:
    input:
        expand('{ref}-sequences.fasta', ref=ref),
        expand('{sample}_{hap}-{ref}-motif_cnt.txt', ref=ref, hap=['hap1', 'hap2'], sample=manifest_df.index)

rule make_paf:
    input:
        fasta = get_fasta,
        ref = ref_path
    output:
        paf = temp('{sample}/{ref}/{sample}_{hap}.paf')
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.24",
    threads: 8
    resources:
        mem = 12,
        hrs = 72
    shell:
        '''
        minimap2 -c -t {threads} -K {resources.mem}G --eqx --cs {minimap_params} --secondary=no --eqx -Y {input.ref} {input.fasta} > {output.paf}
        '''

rule subset_paf:
    input:
        paf = rules.make_paf.output.paf,
        bed = target_region
    output:
        subset = temp('{sample}/{ref}/{sample}_{hap}-subset.paf')
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.29",
    threads: 1
    resources:
        mem = 1,
        hrs = 72
    shell:
        '''
        rustybam liftover --bed {input.bed} {input.paf} > {output.subset}
        '''

rule get_id:
    input:
        subset_paf = rules.subset_paf.output.subset,
    output:
        id = temp('{sample}/{ref}/{sample}_{hap}-id.txt')
    threads: 1
    resources:
        mem = 1,
        hrs = 72
    shell:
        '''
        cut -f1,3,4 {input.subset_paf} | sed 's/\\t/:/; s/\\t/-/' > {output.id}
        '''

rule subset_fasta:
    input:
        fasta = get_fasta,
        id = rules.get_id.output.id
    output:
        subset = '{sample}/{ref}/{sample}_{hap}.fasta'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.14",
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        samtools faidx {input.fasta} $(cat {input.id}) > {output.subset}
        # add sample-hap to header
        sed -i "1 s/$/_{wildcards.sample}-{wildcards.hap}/" {output.subset}
        '''

rule get_motif_count:
    input:
        fasta = rules.subset_fasta.output.subset
    output:
        motif_txt = '{sample}_{hap}-{ref}-motif_cnt.txt',
        mer_cnts_jf = temp('{sample}/{ref}/{hap}-mer_counts.jf')
    params:
        periodicity = len(motif)
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "jellyfish/2.3.0",
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        jellyfish count -m {params.periodicity} -s 100M -t {threads} -C {input.fasta} -o {output.mer_cnts_jf}
        jellyfish query {output.mer_cnts_jf} {motif} > {output.motif_txt}
        '''

rule combine_fasta:
    input:
        fasta = expand('{sample}/{{ref}}/{sample}_{hap}.fasta', hap=['hap1', 'hap2'], sample=manifest_df.index)
    output:
        aggregated = '{ref}-sequences.fasta'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.14",
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        cat {input.fasta} > {output.aggregated}
        '''

rule get_kmer_abundance:
    input:
        agg_fasta = rules.combine_fasta.output.aggregated
    output:
        kmer_abundance = ''
    params:
        periodicity = len(motif)
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "jellyfish/2.3.0",
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        
        '''