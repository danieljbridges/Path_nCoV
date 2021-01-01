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
output_dir = os.path.join("5_GISAID")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
print("  Rampart directory: %s" % rampart_dir)
print("  Artic directory: %s" % artic_dir)
print("  Output directory: %s" % output_dir)
print("Done.")
print("")


EXAMINE CONTENTS
print("Examining directory contents...")
contents_dt = {}
rs = os.listdir(artic_dir)
rs.remove("C02")
for r in os.listdir(artic_dir):
    if os.path.isdir(os.path.join(artic_dir, r)) and r.startswith("C"):
        if r == "C02":
            d = os.path.join(artic_dir, r, "Artic/processed")
        else:
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


# COMPUTE GISAID STATISTICS
print("Computing GISAID statistics...")

# Prepare storage
dts = []

# Iterate over runs
print("  Run  Samples complete")
for r in rs:
    
    # Define run directory
    run_dir = os.path.join(artic_dir, r)
    if os.path.isdir(run_dir) and r.startswith("C"):
        if r == "C02":
            d = os.path.join(run_dir, "Artic/processed")
            slist_fn = os.path.join(rampart_dir, "%s_Artic_SList.txt" % r)
        else:
            d = os.path.join(run_dir, "processed")
            slist_fn = os.path.join(rampart_dir, "%s_SList.txt" % r)
            
        # Load sample-barcode map
        try:
            sample_dt = load_sample_list(slist_fn)
        except:
            print("Failed to load sample list for run %s" % r)


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
            stats_dt["run"] = r
            stats_dt["sample"] = s
            stats_dt["barcode"] = b

            # Calc. most statistics from coverage file
            coverage_df = load_coverage_files(sample_dir)
            stats_dt.update(calc_gisaid_stats(coverage_df))
            
            # Calc. Ns per 100kbp from consensus FASTA
            consensus_path = os.path.join(sample_dir, "%s.consensus.fasta" % s)
            ns_per_100kbp = calc_ns_per_100kbp(consensus_path, verbose=False)
            stats_dt.update({"ns_per_100kbp": ns_per_100kbp})
            
            # Calc. sequencing depth from .fastq
            b_fn = os.path.join(d.replace("processed", "fastq"), "C%d_barcode%.2d.fastq" % (int(r[1:]), b))
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


# WRITE RESULTS
print("Writing results...")
output_fn = "%s_gisaid.csv" % datetime.datetime.now().strftime("%Y-%m-%d")
df.to_csv(os.path.join(output_dir, output_fn), index=False)
print("  To: %s" % os.path.join(output_dir, output_fn))
print("Done.")
print("")


print("-" * 80)
print("Runtime: %s" % str(datetime.timedelta(seconds=time.time() - start_time)))
print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("=" * 80)