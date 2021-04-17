#!/usr/bin/env python
# coding: utf-8
# Merge datasets from sequencing runs
# ----------------------------------------

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
gisaid = pd.read_csv(os.path.join(data_dir, "gisaid.csv"))
#Remove any duplicates from the list
gisaid_df = gisaid[gisaid['Duplicates'] == False]
pango_df = pd.read_csv(os.path.join(data_dir, "lineage_report.csv"))
nextclade_df = pd.read_csv(os.path.join(data_dir, "nextclade.csv"), sep=';')
print("Done.")
print("")

print("Merging gisaid (n=%d), pango (n=%d) and nextclade (n=%d) data together" % (gisaid_df.shape[0], pango_df.shape[0], nextclade_df.shape[0]))
Merge1 = pd.merge(gisaid_df, pango_df, how='inner', left_on='SeqID', right_on='taxon')
Merge2 = pd.merge(Merge1, nextclade_df, how='inner', left_on='SeqID', right_on='seqName')
print("Total remaining records = %d" % (Merge2.shape[0]))
print("Done.")
print("")

print("Removing unwanted column headings")
dropcols = ['PCRPrimers','Type','Duplicates','ExcludeSample','ref_genome_length','pangoLEARN_version','taxon','seqName','qc.mixedSites.mixedSitesThreshold','qc.mixedSites.score','qc.mixedSites.status','qc.mixedSites.totalMixedSites','qc.privateMutations.cutoff','qc.missingData.missingDataThreshold','totalNonACGTNs','nonACGTNs']

Summary_df = Merge2.drop(dropcols, axis=1)
output_fn = "Sequencing_Summary.csv"
Summary_df.to_csv(os.path.join(data_dir, output_fn), index=False)
print("Done.")
print("")
