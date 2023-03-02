#!/usr/bin/env python3

"""
Usage: bed_append.py b1.bed b2.bed
"""

import argparse
import sys
from pybedtools import BedTool
import pandas as pd


def parse_b_beds(fn):
    df = pd.read_csv(fn, sep='\t')

    essential_bed = df.iloc[:, :3].copy()
    df = df.iloc[:,3:]

    def reformat(row):
        return row.index[0:].astype(str) + '=' + row[0:].astype(str)

    essential_bed['info'] = df.apply(lambda x: ';'.join(reformat(x)), axis=1)
    del df

    return essential_bed.sort_values(by=essential_bed.columns[:2].tolist())


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    # parser.add_argument('a', metavar='a.bed',
    #                     help='query bed')
    parser.add_argument('b', metavar='b1.bed b2.bed', nargs='*',
                        help='annotating databases in bed format')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')
    return parser


def main():
    """Get bed"""
    parser = get_parser()
    args = parser.parse_args()

    for b in args.b:
        parse_b_beds(b).to_csv(args.output, sep='\t', index=False, header=True)


if __name__ == "__main__":
    sys.exit(main())
