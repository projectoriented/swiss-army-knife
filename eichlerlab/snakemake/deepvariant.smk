import pandas as pd


configfile: 'config.yaml'
manifest = config.get('manifest', 'manifest.tab')
ref=config.get('ref')

manifest_df = pd.read_csv(manifest, sep='\t', index_col='sample')

chrom_list = ['chr{}'.format(x) for x in list(range(1, 23)) + ['X', 'Y']]

containers={
    'samtools': 'docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1',
    'deepvariant': 'docker://google/deepvariant:1.3.0',
    'bcftools': 'docker://quay.io/biocontainers/bcftools:1.12--h45bccc9_1'
}

wildcard_constraints:
    sample='|'.join(manifest_df.index),
    tech='hifi|ont',
    chr='|'.join(chrom_list)


def get_bam(what_type='bam'):
    def inner(wildcards):
        bam_file_path= manifest_df.at[wildcards.sample, 'bam']
        if what_type == 'bam':
            return bam_file_path
        elif what_type == 'bai':
            return bam_file_path+'.bai'
    return inner

rule all:
    input:
        expand('{sample}/deepvariant/{sample}-{tech}.{ext}', sample=manifest_df.index, tech='hifi', ext=['vcf.gz', 'norm.vcf.gz', 'vcf.gz.tbi', 'norm.vcf.gz.tbi'])

rule split_bam_by_chrom:
    input:
        bam = get_bam(what_type='bam'),
        bai = get_bam(what_type='bai')
    output:
        chr_bam = temp('{sample}/{chr}-{tech}.bam'),
        chr_bam_bai = temp('{sample}/{chr}-{tech}.bam.bai'),
    threads: 4
    resources:
        mem = lambda wildcards, attempt: attempt * 1,
        hrs=12,
    container: containers['samtools']
    shell:
        '''
        samtools view -@ {threads} -h -b {input.bam} -o {output.chr_bam} {wildcards.chr}
        samtools index -@ {threads} {output.chr_bam} -o {output.chr_bam_bai}
        '''

rule deepvariant_by_chrom:
    input:
        ref = ref,
        bam = rules.split_bam_by_chrom.output.chr_bam,
        bai= rules.split_bam_by_chrom.output.chr_bam_bai
    output:
        vcf_gz = temp('{sample}/deepvariant/{sample}-{tech}-{chr}.vcf.gz'),
    threads: 8
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (2 * threads),
        hrs=12,
    container: containers['deepvariant']
    benchmark: 'benchmarks/{sample}_{tech}-deepvariant_benchmark-{chr}.txt'
    shell:
        '''
        /opt/deepvariant/bin/run_deepvariant \
            --model_type="PACBIO" \
            --ref={input.ref} \
            --reads={input.bam} \
            --output_vcf={output.vcf_gz} \
            --call_variants_extra_args="use_openvino=true" \
            --num_shards={threads} \
            --dry_run=false
        '''

rule vcf_norm:
    input:
        vcf = rules.deepvariant_by_chrom.output.vcf_gz,
        ref = ref
    output:
        norm_vcf = temp('{sample}/deepvariant/{sample}-{tech}-{chr}.norm.vcf.gz')
    threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (threads),
        hrs=72
    container: containers['bcftools']
    shell:
        '''
        bcftools norm -f {input.ref} {input.vcf} -Oz -o {output.norm_vcf}
        '''

rule gather_vcfs:
    input:
        vcfs = expand('{{sample}}/deepvariant/{{sample}}-{{tech}}-{chr}.vcf.gz', chr=chrom_list),
        norm_vcfs = expand('{{sample}}/deepvariant/{{sample}}-{{tech}}-{chr}.norm.vcf.gz', chr=chrom_list),
    output:
        merged_vcf_gz = '{sample}/deepvariant/{sample}-{tech}.vcf.gz',
        merged_norm_vcf_gz = '{sample}/deepvariant/{sample}-{tech}.norm.vcf.gz',
    threads: 1
    resources:
        mem = lambda wildcards, attempt: attempt * 2,
        hrs=72
    container: containers['bcftools']
    shell:
        '''
        bcftools concat {input.vcfs} -Oz -o {output.merged_vcf_gz}
        bcftools concat {input.norm_vcfs} -Oz -o {output.merged_norm_vcf_gz}
        '''

rule tabix:
    input:
        vcf = rules.gather_vcfs.output.merged_vcf_gz,
        norm_vcf = rules.gather_vcfs.output.merged_norm_vcf_gz,
    output:
        vcf_tbi = '{sample}/deepvariant/{sample}-{tech}.vcf.gz.tbi',
        vcf_norm_tbi = '{sample}/deepvariant/{sample}-{tech}.norm.vcf.gz.tbi',
    threads: 1
    resources:
        mem = lambda wildcards, attempt: attempt * 2,
        hrs=12,
    container: containers['samtools']
    shell:
        '''
        tabix -p vcf {input.vcf}
        tabix -p vcf {input.norm_vcf}
        '''