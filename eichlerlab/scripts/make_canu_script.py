#!/usr/bin/env python3

# Usage: make_canu_script.py --sample 14455_p1 --trio
# NOTE if you're using the trio flag, and you have a different Illumina directory- feel free to change it.

from glob import glob, iglob
import re
import pandas as pd
import sys
import argparse
import os

ILLUMINA_READS_DIR = '/path/to/Illumina/WGS'
LRA_DIR = '/path/to/clinical'


def find_illumina_fofn(s, prefix):
    fp = [x for x in iglob(f'{prefix}/**/*.fofn', recursive=True) if s in x]
    if fp:
        return fp[0]
    else:
        sys.exit(f'No Illumina fofn for {s}')


def make_ont_fofn(sample, fn, prefix):
    pattern = r'(.*fastq\/)(.*)(\/guppy\/)(.*)\/(sup-prom)\/(.*fastq_pass\.fastq\.gz)'
    a = []
    file_list = glob(
        f'{prefix}/{sample}/raw_data/nanopore/STD/**/*fastq_pass.fastq.gz',
        recursive=True)
    for f in file_list:
        match = re.findall(pattern, f)
        if match:
            d = {}
            match = match[0]
            d['sample'] = sample
            d['flow_cell'] = match[1]
            d['version'] = match[3]
            d['model'] = match[4]
            d['path'] = f
            a.append(d)

    fofn_df = pd.DataFrame(a)
    fofn_df.drop_duplicates(subset=['sample', 'flow_cell'], keep='last', inplace=True)
    return fofn_df['path'].to_csv(fn, header=False, index=False)


def make_script(sample, illumina_prefix, trio=False):
    family_id = re.findall(r'[a-zA-Z0-9]+', sample)[0]

    mat_fofn = find_illumina_fofn(s=f'{family_id}_mo', prefix=illumina_prefix)
    pat_fofn = find_illumina_fofn(s=f'{family_id}_fa', prefix=illumina_prefix)

    if trio:
        canu_script = (
            "#!/bin/bash\n"
            "module load canu/2.1.1\n"
            "canu -p {s} -haplotype genomeSize=3.1g useGrid=true maxInputCoverage=150 gridEngineResourceOption='-pe serial THREADS -l mfree=MEMORY' -haplotype1 $(cat {pat_fp}) -haplotype2 $(cat {mat_fp}) -nanopore $(cat {ont_fp}) maxThreads=16 maxMemory=128G cnsMemory=128G -d {dn}"
        ).format(s=sample, pat_fp=pat_fofn, mat_fp=mat_fofn,
                 dn=f'{sample}-ont', ont_fp=f'{sample}-ont-guppy.fofn')
    else:
        canu_script = (
            "#!/bin/bash\n"
            "module load canu/2.1.1\n"
            "canu -p {s} -d {dn} genomeSize=3.1g maxInputCoverage=150 gridEngineResourceOption='-pe serial THREADS -l mfree=MEMORY' -nanopore {fp} maxThreads=16 maxMemory=128G cnsMemory=128G"
        ).format(s=sample, dn=f'{sample}-ont', fp=f'{sample}-ont-guppy.fofn')

    return canu_script


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', '-s', type=str, required=True,
                        help='Supply sample(s) e.g. 14455_p1 14455_mo')
    parser.add_argument('--trio', required=False, help='Get script for trio. (default: False)',
                        action='store_true', default=False)
    parser.add_argument('--lra', '-l', type=str, required=True,
                        help='Parent directory of the long read archive')
    parser.add_argument('--illumina', '-m', type=str, required=True,
                        help='Parent directory of Illumina fastqz')
    return parser


def main():
    """Make canu script"""
    parser = get_parser()
    args = parser.parse_args()

    sample = args.sample

    # Fixate prefix.
    prefix_lra = args.lra if args.lra else LRA_DIR
    prefix_illumina = args.illumina if args.illumina else ILLUMINA_READS_DIR

    canu_script = make_script(sample, trio=args.trio, illumina_prefix=prefix_illumina)
    canu_script_fn = f'run-canu_{sample}.sh'
    if args.trio:
        canu_script_fn = f'run-canu_{sample}-trio.sh'

    # Write ont-fofn.
    ont_fofn = f'{sample}-ont-guppy.fofn'
    make_ont_fofn(sample, ont_fofn, prefix=prefix_lra)

    # Write canu script
    with open(canu_script_fn, 'w+') as f:
        f.writelines(canu_script)

    os.chmod(canu_script_fn, 0o775)


if __name__ == "__main__":
    sys.exit(main())
