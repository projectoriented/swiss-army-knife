#!/usr/bin/env python3

"""
Usage: bed_append.py -a b1.bed -b b2.bed b3.bed --bn b2 b3
"""

import argparse
import sys
from pybedtools import BedTool
import pandas as pd
import re
from itertools import chain


def parse_b_beds(fn):
    with open(fn, 'r') as f:
        match = re.findall(r'(chr[0-9]+)\t', f.readline())
    header = 0 if not match else None
    df = pd.read_csv(fn, sep='\t', header=header)

    essential_bed = df.iloc[:, :3].copy()
    df = df.iloc[:, 3:]

    def reformat(row):
        if header == 0:
            return ';'.join(row.index[0:].astype(str) + '=' + row[0:].astype(str))
        else:
            return ((len(row) - 1) * "{};" + "{}").format(*row.tolist())

    essential_bed['info'] = df.apply(lambda x: reformat(x), axis=1)
    del df

    # Sort
    essential_bed = essential_bed.sort_values(by=essential_bed.columns[:2].tolist())

    return essential_bed


def add_to_query_bed(a_df, b_bed, b_bed_name):
    a_headers = a_df.columns.str.upper().tolist()

    a_bed_obj = BedTool.from_dataframe(a_df)

    b_bed = b_bed.sort()
    intersected = a_bed_obj.intersect(b_bed, loj=True, wa=True)

    ab_headers = a_headers + list(chain(*map(lambda x: ['chrom' + x, 'start' + x, 'end' + x, x], [b_bed_name])))

    intersected_df = intersected.to_dataframe(names=ab_headers)

    intersected_df = intersected_df[a_headers + [b_bed_name]].sort_values(intersected_df.columns[:2].tolist())
    intersected_df.drop_duplicates(subset=intersected_df.columns[:2], inplace=True)
    intersected_df.set_index(a_headers, inplace=True)

    return intersected_df


def bedtools_intersect(a_df, b_bed):
    a_bed_obj = BedTool.from_dataframe(a_df)
    b_bed = b_bed.sort()
    intersected = a_bed_obj.intersect(b_bed, loj=True, wa=True)
    intersected_df = intersected.to_dataframe()
    result_column = intersected_df.iloc[:, -1:]
    return result_column


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', metavar='a.bed',
                        help='query bed')
    parser.add_argument('-b', metavar='b1.bed b2.bed', nargs='*',
                        help='annotating databases in bed format')
    parser.add_argument('--bn', metavar='b1 b2', nargs='*',
                        help='b bed name delimited by space and in the same order of b bed files')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')
    return parser


def main():
    """Get bed"""
    parser = get_parser()
    args = parser.parse_args()

    df_a = pd.read_csv(args.a, sep='\t', header=0, low_memory=False)
    essential_df = df_a.iloc[:, :3].copy()
    df_a = df_a.iloc[:, 3:]

    bed_names = args.bn

    final_df = pd.DataFrame()
    for idx, b in enumerate(args.b):
        print(f'parsing {b}')
        bed_object = BedTool.from_dataframe(parse_b_beds(b))
        intersected = add_to_query_bed(a_df=essential_df, b_bed=bed_object, b_bed_name=bed_names[idx])
        final_df = pd.concat([final_df, intersected], axis=1)
        print(f'intersected {b}')

    final_df.reset_index(inplace=True, drop=True)
    final_df = pd.concat([essential_df, df_a, final_df], axis=1)
    final_df.to_csv(args.output, sep='\t', header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
