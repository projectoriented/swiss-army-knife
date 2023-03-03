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
            return ((len(row)-1)*"{};" + "{}").format(*row.tolist())

    essential_bed['info'] = df.apply(lambda x: reformat(x), axis=1)
    del df

    # Sort
    essential_bed = essential_bed.sort_values(by=essential_bed.columns[:2].tolist())

    return essential_bed


def add_to_query_bed(a_df, b_bedtool_objects: list, b_bed_names: list):
    a_headers = a_df.columns.str.upper().tolist()

    a_bed_obj = BedTool.from_dataframe(a_df)
    # b_bed_obj = BedTool.from_dataframe(df)

    intersected = a_bed_obj.intersect(b_bedtool_objects, wa=True, wb=True)

    ab_headers = a_headers + list(chain(*map(lambda x: ['chrom' + x, 'start' + x, 'end' + x, x], b_bed_names)))

    intersected_df = intersected.to_dataframe(names=ab_headers)
    # intersected_df = pd.DataFrame.from_records([f[0:len(a_headers)] + f[-1:] for f in intersected], columns=a_headers+[b_bed_name])

    return intersected_df[a_headers + b_bed_names].sort_values(intersected_df.columns[:2].tolist())


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
    bed_names = args.bn

    b_dict = {}
    for idx, b in enumerate(args.b):
        b_bt_object = BedTool.from_dataframe(parse_b_beds(b))
        b_dict[bed_names[idx]] = b_bt_object

    b_df_collection = list(b_dict.values())
    final_df = add_to_query_bed(a_df=df_a, b_bedtool_objects=b_df_collection, b_bed_names=bed_names)
    final_df.to_csv(args.output, sep='\t', header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
