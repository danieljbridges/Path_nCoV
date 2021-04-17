#!/usr/bin/env python3
# Compute statistics for GISAID submission
# ----------------------------------------
#
# Example usage:
#   python run_gisaid-statistics.py -d data/WGS
#
# JHendry, 2021/01/01


import os
import sys
import subprocess
import time
import datetime
import getopt
import pandas as pd
import numpy as np
from lib.gisaid import *


print("=" * 80)
print("Compute GISAID statistics")
print("-" * 80)
print("Command: %s" % " ".join(sys.argv))
print("Run on host: %s" % os.uname().nodename)
print("Operating system: %s" % os.uname().sysname)
print("Machine: %s" % os.uname().machine)
print("Started at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("=" * 80)
start_time = time.time()


# PARSE CLI INPUT
print("Parsing command line inputs...")
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


# PREPARE DIRECTORIES
print("Preparing directories...")
rampart_dir = os.path.join(data_dir, "2_SampleList_and_Rampart")
artic_dir = os.path.join(data_dir, "3_Artic_Output")
output_dir = os.path.join(data_dir, "5_GISAID")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
print("  Rampart directory: %s" % rampart_dir)
print("  Artic directory: %s" % artic_dir)
print("  Output directory: %s" % output_dir)
print("Done.")
print("")


# EXAMINE CONTENTS
print("Examining directory contents...")
contents_dt = {}
rs = os.listdir(artic_dir)
for r in os.listdir(artic_dir):
    if os.path.isdir(os.path.join(artic_dir, r)) and r.startswith("C"):
        d = os.path.join(artic_dir, r, "processed")
        n_samples = sum([1 for s in os.listdir(d) 
                         if os.path.isdir(os.path.join(d, s))])
        contents_dt[r] = n_samples
print("  Run\tNo. samples")
for d, n in contents_dt.items():
    print("  %s\t%d" % (d, n))
print("  Total runs: %d" % len(contents_dt.keys()))
print("  Total samples: %d" % sum(contents_dt.values()))
print("Done.")
print("")


# LOAD SAMPLE METADATA
print("Checking sample processing details...")
sample_df = pd.read_csv(os.path.join(rampart_dir, "Samples_Sequenced.csv"))
print("  Total runs: %d" % sample_df["SeqRun"].unique().shape[0])
print("  Total samples: %d" % sample_df["SeqID"].unique().shape[0])
print("  Unique samples: %d" % sample_df["SampleID"].unique().shape[0])
print("  Barcodes used: %d" % sample_df["SeqBarcode"].unique().shape[0])
# Map from Unique ID to barcode ID
sample_dt = { row["SeqID"]: row["SeqBarcode"] for _, row in sample_df.iterrows() }
print("Done.")
print("")


# COMPUTE SEQUENCING STATISTICS
print("Computing sequencing statistics...")

# Prepare storage
dts = []

# Iterate over runs
print("  Run  Samples/Complete")
for r in rs:
    
    # Define run directory
    run_dir = os.path.join(artic_dir, r)
    if os.path.isdir(run_dir) and r.startswith("C"):
        d = os.path.join(run_dir, "processed")

        # Iterate over samples
        n_total = len(os.listdir(d))
        for i, s in enumerate(os.listdir(d)):
            
            # Print progress
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write("  %s %d/%d" % (r, i+1, n_total))
            
            # Define sample directory
            sample_dir = os.path.join(d, s)
            
            # Define sample barcode
            b = sample_dt[s]
            
            # Identifiers
            stats_dt = {}
            stats_dt["SeqRun"] = r
            stats_dt["SeqID"] = s
            stats_dt["SeqBarcode"] = b

            # Calc. most statistics from coverage file
            coverage_df = load_coverage_files(sample_dir)
            stats_dt.update(calc_gisaid_stats(coverage_df))
            
            # Calc. Ns per 100kbp from consensus FASTA
            consensus_path = os.path.join(sample_dir, "%s.consensus.fasta" % s)
            ns_per_100kbp = calc_ns_per_100kbp(consensus_path, verbose=False)
            stats_dt.update({"ns_per_100kbp": ns_per_100kbp})
            
            # Calc. sequencing depth from .fastq
            b_fn = os.path.join(d.replace("processed", "fastq"), "C%02d_barcode%02d.fastq" % (int(r[1:]), b))
            try:
                sequencing_depth_avg_fastq = calc_avg_seq_depth(b_fn, genome_length=stats_dt["ref_genome_length"])
                stats_dt.update({"sequencing_depth_avg_fastq": sequencing_depth_avg_fastq})
            except:
                stats_dt.update({"sequencing_depth_avg_fastq": 0})
            
            # Store
            dts.append(stats_dt)
        
        sys.stdout.flush()
        sys.stdout.write("\n")
        
# Create data frame
df = pd.DataFrame(dts)
print("Done.")
print("")


# MERGE WITH SAMPLE LIST
print("Merging results with sample list...")
print("  No. samples...")
print("    ...in sample list: %d" % sample_df.shape[0])
print("    ...with consensus sequence: %d" % df.shape[0])
# Merge
merged_df = pd.merge(left=sample_df,
                     right=df,
                     left_on=["SeqRun", "SeqBarcode", "SeqID"],
                     right_on=["SeqRun", "SeqBarcode", "SeqID"])
print("    ...after merging: %d" % merged_df.shape[0])
print("Done.")
print("")

# Highlight duplicates
print("Highlighting highest depth for duplicate samples...")
grps = merged_df.groupby("SampleID")

l_dfs = []
print("  {:<8}  {:<8}  {:<10}  {:<4}".format("Sample", "No. dup.", "Keep", "Depth"))
for n, sdf in grps:
    n_dup = sdf.shape[0]
    keep = sdf.sort_values("sequencing_depth_avg", ascending=False)
    values = sdf.sort_values("sequencing_depth_avg", ascending=False).iloc[0]
    keep["Duplicates"] = [False]+[True]*(n_dup-1)
    if n_dup > 1 :
        print("  {:<8}  {:<8} {:<25} {:<4.1f}".format(n, n_dup, values['SeqID'], values['sequencing_depth_avg']))
    l_dfs.append(keep)
retained_df = pd.concat(l_dfs, 0)
print("  Samples remaining: %d" % retained_df.shape[0])
print("Done.")
print("")

# WRITE RESULTS
print("Writing results...")
output_fn = "gisaid.csv"
retained_df.to_csv(os.path.join(output_dir, output_fn), index=False)
print("  To: %s" % os.path.join(output_dir, output_fn))
print("Done.")
print("")


print("-" * 80)
print("Runtime: %s" % str(datetime.timedelta(seconds=time.time() - start_time)))
print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("=" * 80)
