#!/usr/bin/env python3

"""
Usage: ./filter_data_tbl.py variants_asd_families_sv_insdel.tsv.gz
"""

import sys
import pandas as pd


def insert_str(string, str_to_insert, index):
    return string[:index] + str_to_insert + string[index:]


original_data_tbl = sys.argv[1]
insert_pos = original_data_tbl.find('.')
out_fn = insert_str(original_data_tbl, '_filt', insert_pos)

df = pd.read_csv(original_data_tbl, sep='\t', header=0)

# Write out
df[(~df['PBSV'].isna()) & (df['VALID'] == 'VALID')].to_csv(out_fn, sep='\t', header=True, index=False)
print(f'file written to: {out_fn}')
