#!/usr/bin/env python
# coding: utf-8
# Merge datasets from sequencing runs
# ----------------------------------------
# Example usage:
#   python run_gisaid-statistics.py -d data/WGS
#
# JHendry, 2021/01/01

import os
import pandas as pd
import numpy as np
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], ":d:")
except getopt.GetoptError:
    print("  Error parsing options.")
    sys.exit(2)
for opt, value in opts:
    if opt == "-d":
        data_dir = value
        print("  Data directory: %s" % data_dir)
    else:
        print("  Parameter %s not recognized." % opt)
        sys.exit(2)
print("Done.")
print("")

print("Reading in gisaid, pangolin and nextclade arguments.")
gisaid_df = pd.read_csv(os.path.join(data_dir, "gisaid.csv"))
pango_df = pd.read_csv(os.path.join(data_dir, "lineage_report.csv"))
nextclade_df = pd.read_csv(os.path.join(data_dir, "nextclade.csv"), sep=';')
print("Done.")
print("")

print("Merging datasets together")
Merge1 = pd.merge(gisaid_df, pango_df, how='inner', left_on='sample', right_on='taxon')
Merge2 = pd.merge(Merge1, nextclade_df, how='inner', left_on='sample', right_on='seqName')
print("Done.")
print("")

print("Removing unwanted column headings")
dropcols = ['Primers','Type','ref_genome_length','pangoLEARN_version','taxon','seqName']
Summary_df = Merge2.drop(dropcols, axis=1)
output_fn = "Summary.csv"
Summary_df.to_csv(os.path.join(data_dir, output_fn), index=False)
print("Done.")
print("")
