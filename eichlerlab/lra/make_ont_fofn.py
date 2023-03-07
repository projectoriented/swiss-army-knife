#!/usr/bin/env python3

# Usage: ./make_ont_fofn.py --sample 14455_p1 --proj_dir path/to/clinical --output 14455_p1.fofn
# This script will get you all the fastqs that belong to a major version of guppy for STD library
# future implementations will include: --filter_string 'lib=STD;model=sup-prom;version=6.1.2'
# Slack/email me if you experience any problems or want additional features, wumei@uw.edu

from glob import iglob
import re
import pandas as pd
import sys
import argparse
import os


def make_ont_fofn(sample, fn, prefix):
    pattern = r'(.*fastq\/)(.*)(\/guppy\/)(.*)\/(sup-prom)\/(.*fastq_pass\.fastq\.gz)'
    a = []
    search_path = os.path.join(prefix, sample, 'raw_data/nanopore/STD/**/*fastq_pass.fastq.gz')

    file_list = iglob(search_path, recursive=True)
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


def main():
    """Run script to output ont fofn"""
    parser = get_parser()
    args = parser.parse_args()

    make_ont_fofn(sample=args.sample, fn=args.output, prefix=args.proj_dir)


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument('--sample', type=str, required=True,
                        help='e.g. 14455_p1')
    parser.add_argument('--proj_dir', type=str, required=True,
                        help='Absolute path for project directory where fast5s are found.')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')

    return parser


if __name__ == "__main__":
    sys.exit(main())
