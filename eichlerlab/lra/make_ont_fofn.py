#!/usr/bin/env python3

# Usage: ./make_ont_fofn.py --sample 14455_p1 --proj_dir absolute/path/to/clinical --output 14455_p1.fofn --filter_string 'lib=STD;model=sup-prom;bc=guppy;ver=6'
# This script will get you all the fastqs that belong to guppy ver 6.x.x for STD library
# Slack/email me if you experience any problems or want additional features, wumei@uw.edu

from glob import iglob
import re
import pandas as pd
import sys
import argparse
import os


def make_ont_fofn(sample, fn, prefix, fltr_str):
    fltr_dict = {
        'bc': 'guppy',
        'model': 'sup-prom',
        'ver': '6',
        'lib': 'STD'
    }
    if fltr_str:
        alt_fltr_list = fltr_str.split(';')
        alt_fltr_dict = dict([x.split('=') for x in alt_fltr_list])
        fltr_dict.update(alt_fltr_dict)

    basecaller = fltr_dict['bc']
    model = fltr_dict['model']
    version = fltr_dict['ver']
    library = fltr_dict['lib']

    regex = fr'(.*nanopore)/(?P<lib>{library})/(.*fastq)/(?P<runid>.*)/(?P<bc>{basecaller})/(?P<ver>{version}.*)/(?P<model>{model})/(?P<fn>.*fastq_pass\.fastq\.gz)'
    search_pattern = re.compile(regex)
    a = []
    search_path = os.path.join(prefix, sample, 'raw_data/nanopore/*/**/*fastq_pass.fastq.gz')

    file_list = iglob(search_path, recursive=True)
    for f in file_list:
        match = search_pattern.match(f)
        if match:
            match_dict = match.groupdict()
            match_dict['fpath'] = f
            a.append(match_dict)

    fofn_df = pd.DataFrame(a)
    if fofn_df.empty:
        sys.exit(f'No fastq.gz matched the filter- please try again :)')
    return fofn_df['fpath'].to_csv(fn, header=False, index=False)


def main():
    """Run script to output ont fofn"""
    parser = get_parser()
    args = parser.parse_args()

    make_ont_fofn(sample=args.sample, fn=args.output, prefix=args.proj_dir, fltr_str=args.filter_string)


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
    parser.add_argument('--filter_string', type=str, required=False,
                        help='Absolute path for project directory where fast5s are found.')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')

    return parser


if __name__ == "__main__":
    sys.exit(main())
