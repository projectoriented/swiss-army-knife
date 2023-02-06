#!/usr/bin/env python3
"""
Get post-assembly stats for a sample.
Usage: ./post_assembly_stats.py --sample 14455_p1
"""

import os
import sys
import pandas as pd
import argparse
import re
import numpy as np


def get_merqury_stats_df(prefix, s):
    # sample name workaround right now
    sample = f'{s}_HIFI'

    def hapmer_resolved():
        """
        Check if merqury was run with hapmer information.
        :return: parent directory path
        """
        trio_dir = os.path.join(prefix, sample, 'trio')
        standard = os.path.split(trio_dir)[0]
        if os.path.exists(trio_dir):
            return trio_dir
        elif os.path.exists(standard):
            return standard
        else:
            sys.exit(f'No merqury files found for {sample}')

    working_path = hapmer_resolved()

    def get_fn(which_fn='qv'):
        """
        Using the parent directory, find the target files.
        :param which_fn: qv or completeness
        :return: the full path of the target file
        """
        if which_fn == 'qv':
            if 'trio' in working_path:
                return os.path.join(working_path, f'{sample}_trio.qv')
            else:
                return os.path.join(working_path, f'{sample}.qv')
        elif which_fn == 'completeness':
            return os.path.join(working_path, 'completeness.stats')

    df_qv = pd.read_csv(get_fn(which_fn='qv'), sep='\t', header=None, nrows=2,
                        names=['haplotype', 'assembly_kmers', 'both_kmers', 'qv', 'error_rate'])
    df_qv = df_qv.astype({'qv': int})
    qv = df_qv['qv'].astype(str).str.cat(sep='|')  # condense to reflect hap1|hap2
    del df_qv

    df_c = pd.read_csv(get_fn(which_fn='completeness'), sep='\t', header=None,
                       names=['haplotype', 'both_kmers', 'assembly_kmers', 'read_kmers', 'completeness'])
    df_c = df_c.astype({'completeness': int})
    if 'trio' in working_path:
        df_c = df_c.iloc[7:9]
    else:
        df_c = df_c.iloc[0:2]
    c = df_c['completeness'].astype(str).str.cat(sep='|')  # condense to reflect hap1|hap2
    del df_c

    return pd.DataFrame({'completeness': c, 'qv': qv}, index=[s])


def get_n50(prefix, sample):
    phase = 'dip'
    if ('mo' in sample) or ('fa' in sample):
        phase = 'bp'

    parent_dir = os.path.join(prefix, f'{sample}/assemblies/hifiasm/')
    try:
        target_fn = [target for target in os.listdir(parent_dir) if phase in target][0]
        fp = os.path.join(parent_dir, target_fn)
    except IndexError:
        sys.exit(f'No n50 file exists for {sample}')

    def parse_fn(_):
        assembler_ver = re.findall(r'\d+\.\d+\.\d+', target_fn)[0]
        return phase, assembler_ver

    n50_dict = {}
    with open(fp) as fp:
        for idx, line in enumerate(fp.readlines()):
            sp_l = "".join(line.split())
            if 'hap1' in sp_l:
                hap1 = sp_l
                n50_dict[hap1] = {}
            elif 'hap2' in sp_l:
                hap2 = sp_l
                n50_dict[hap2] = {}

            entries = sp_l.split(":")
            if (idx < 6) and ('hap1' not in sp_l):
                n50_dict[hap1][entries[0]] = entries[1]
            elif (6 < idx < 13) and ('hap2' not in sp_l):
                n50_dict[hap2][entries[0]] = entries[1]

    phase_n_ver = parse_fn(fp)
    d_keys = list(n50_dict.keys())
    n50 = "|".join([n50_dict[d_keys[0]].get('N50(Mbp)'), n50_dict[d_keys[1]].get('N50(Mbp)')])

    return pd.DataFrame({'N50(Mbp)': n50, 'phase': phase_n_ver[0], 'assembly_ver': phase_n_ver[1]}, index=[sample])


def get_coverage_x(prefix, sample):
    genome_gbp = 3 * (10 ** 9)

    def find_fofn():
        fp = os.path.join(prefix, sample, 'raw_data/PacBio_HiFi/fofn/ccs/fastq.fofn')
        if os.path.exists(fp):
            return fp
        else:
            sys.exit(f'No fofn exists for {sample}')

    fofn_fp = find_fofn()

    fofn_df = pd.read_csv(fofn_fp, sep='\t', header=None, index_col=0)
    df = pd.concat([pd.read_csv(file + '.fai', sep='\t', header=None, usecols=[0, 1]) for file in fofn_df.index])
    len_list = pd.Series(df[1].copy())

    # len_list.sort_values(ascending=False, inplace=True)
    coverage = np.round(np.sum(len_list) / genome_gbp, 2)

    return pd.DataFrame({'Coverage X': coverage}, index=[sample])


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument('--sample', '-s', type=str, required=True,
                        help='Supply sample(s) e.g. 14455_p1 14455_mo')
    parser.add_argument('--header', required=False, help='Get column names for table. (default: False)',
                        action='store_true', default=False)
    parser.add_argument('--merqury', '-m', type=str, required=True,
                        help='Parent directory of Merqury output')
    parser.add_argument('--lra', '-l', type=str, required=True,
                        help='Parent directory of the long read archive')
    parser.add_argument('--assembly', '-a', type=str, required=True,
                        help='Parent directory of hifiasm output')
    parser.add_argument('--output', '-o', type=str, required=False,
                        help='Output tsv to write to',
                        default='/dev/stdout')
    return parser


def main():
    """Run script to get comprehensive post-assembly statistics on a sample"""
    parser = get_parser()
    args = parser.parse_args()

    mq_df = get_merqury_stats_df(args.merqury, args.sample)
    n50_df = get_n50(args.assembly, args.sample)
    cov_df = get_coverage_x(args.lra, args.sample)

    df = pd.concat([mq_df, n50_df, cov_df], axis=1)
    df.index.rename('sample', inplace=True)

    df.to_csv(args.output, sep='\t', header=args.header)


if __name__ == "__main__":
    sys.exit(main())
