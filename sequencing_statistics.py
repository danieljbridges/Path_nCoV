#!/usr/bin/env python
# coding: utf-8
# ----------------------------------------
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
gisaid_dir = os.path.join(data_dir, "5_GISAID")
if not os.path.isdir(gisaid_dir):
    os.makedirs(gisaid_dir)
print("  Rampart directory: %s" % rampart_dir)
print("  Artic directory: %s" % artic_dir)
print("  GISAID directory: %s" % gisaid_dir)
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

# Prepare storage
dts = []

# LOAD SAMPLE METADATA
print("Loading sample list details from Samples_Sequenced.csv...")
samples_df = pd.read_csv(os.path.join(rampart_dir, "Samples_Sequenced.csv"))
print("    %s PCR reactions have been performed consisting of..." % samples_df.shape[0])
print("        %s Controls" % samples_df[samples_df.Type=='Control'].shape[0])
print("        %s Samples, of which...." % samples_df[samples_df.Type=='Sample'].shape[0])
print("            %s are unique" % samples_df[samples_df.Type=='Sample']["SampleID"].unique().shape[0])

print("Removing samples that have not been sequenced")
seqsamples_df = samples_df[samples_df.SeqRun.notnull()]
print("    %s sequencing reactions have been performed consisting of..." % seqsamples_df.shape[0])
print("        %s Controls" % seqsamples_df[seqsamples_df.Type=='Control'].shape[0])
print("        %s Samples, of which...." % seqsamples_df[seqsamples_df.Type=='Sample'].shape[0])
print("            %s are unique" % seqsamples_df[seqsamples_df.Type=='Sample']["SampleID"].unique().shape[0])
# Map from Unique ID to barcode ID
sample_dt = { row["SeqID"]: row["SeqBarcode"] for _, row in seqsamples_df.iterrows() }
print("Done.")
print("")

# COMPUTE SEQUENCING STATISTICS
print("Computing sequencing statistics...")

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
gisaid_df = pd.DataFrame(dts)
print("Done.")

#Check that all the sequence IDs are unique
if len(gisaid_df["SeqID"]) == len(set(gisaid_df["SeqID"])):
    print("    All sequence IDs are unique")
else:
    print("    Duplicate sequence IDs identified")
    l_func = lambda x, y: list((set(x)- set(y))) + list((set(y)- set(x))) 
    non_match = l_func(gisaid_df["SeqID"], gisaid_df["SeqID"].unique()) 
    print("    Non-match elements: ", non_match)
    print("Exiting script")
    exit()
print("")
print("Statistics for a total of %s samples have been calculated" % gisaid_df.shape[0])
gisaid_fn = "intermediates/gisaid.csv"
gisaid_df.to_csv(os.path.join(gisaid_dir, gisaid_fn), index=False)
print("    Data output to %s" % gisaid_fn)
print("Done.")
print("")

#Find the other stats files for merging
print("-" * 80)
print("Combining with other data outputs")
if os.path.isfile(os.path.join(gisaid_dir,"intermediates/lineage_report.csv")) and os.path.isfile(os.path.join(gisaid_dir,"intermediates/nextclade.csv")):
    print ("  PANGO and nextclade files found. Reading in data")
    pango_df = pd.read_csv(os.path.join(gisaid_dir, "intermediates/lineage_report.csv"))
    nextclade_df = pd.read_csv(os.path.join(gisaid_dir, "intermediates/nextclade.csv"), sep=';')
elif os.path.isfile(os.path.join(gisaid_dir,"intermediates/lineage_report.csv")):
    print ("  Only PANGO file found. Exiting script")
    exit()
elif os.path.isfile(os.path.join(gisaid_dir,"intermediates/nextclade.csv")):
    print ("  Only nextclade file found. Exiting script")
    exit()
else:
    print ("  PANGO and nextclade files not found. Exiting script")
    exit()
    
#Check that all sequences are represented. If not this suggests that some of the fasta concatenation has failed or that the wrong sequence name has been used
l_func = lambda x, y: list((set(x)- set(y))) + list((set(y)- set(x))) 
non_match_gp = l_func(gisaid_df["SeqID"], pango_df["taxon"]) 
non_match_gn = l_func(gisaid_df["SeqID"], nextclade_df["seqName"]) 
if len(non_match_gp) > 0 :
    print("  A total of %s sequence ids do not match between gisaid and PANGO" % len(non_match_gp))
    exit()
elif len(non_match_gn) > 0 :
    print("  A total of %s sequence ids do not match between gisaid and nextclade" % len(non_match_gn))
    exit()
else:
    print("  All sequences match across the three files")

print("  Merging gisaid (n=%d), pango (n=%d) and nextclade (n=%d) data together" % (gisaid_df.shape[0], pango_df.shape[0], nextclade_df.shape[0]))
Merge1 = pd.merge(gisaid_df, pango_df, how='inner', left_on='SeqID', right_on='taxon')
Merge2 = pd.merge(Merge1, nextclade_df, how='inner', left_on='SeqID', right_on='seqName')
print("  Total remaining records = %d" % (Merge2.shape[0]))
print("  Removing unwanted column headings from final summary file of samples sequenced")
gisaidcols = ['ref_genome_length']
pangocols = (['pangoLEARN_version','taxon'])
nextcladecols = ['seqName','qc.mixedSites.mixedSitesThreshold','qc.mixedSites.score','qc.mixedSites.status','qc.mixedSites.totalMixedSites','qc.privateMutations.cutoff','qc.missingData.missingDataThreshold','totalNonACGTNs','nonACGTNs']
dropcols = gisaidcols + pangocols + nextcladecols
sequenced_df = Merge2.drop(dropcols, axis=1)
print("Done.")
print("")

# MERGE WITH SAMPLE LIST
print("Merging sequenced samples with sample list...")
print("  No. samples...")
print("    ...in sample list: %d" % seqsamples_df.shape[0])
print("    ...with consensus sequence: %d" % gisaid_df.shape[0])
merged_df = pd.merge(left=seqsamples_df,
                     right=sequenced_df,
                     left_on=["SeqRun", "SeqBarcode", "SeqID"],
                     right_on=["SeqRun", "SeqBarcode", "SeqID"],
                    how='outer')
print("    ...after merging: %d" % merged_df.shape[0])
print("Done.")
print("")

qc_depth = 50
qc_breadth = 50

print("Removing controls, unsequenced samples, those marked for exclusion or failing qc parameters")
keepers_df = merged_df.query("sequencing_depth_avg >= @qc_depth" +
                             "& coverage_breadth >= @qc_breadth" +
                             "& ExcludeSample != 'Y'" +
                             "& Type == 'Sample'" ,
                             inplace=False)
print("    Samples remaining: %d" % keepers_df.shape[0])
print("Done")
print("")

#Highlight duplicates
print("Highlighting highest depth for duplicate samples...")

l_dfs = []
print("  {:<8}  {:<8}  {:<10}  {:<4}".format("Sample", "No. dup.", "Keep", "Depth"))


for n, sdf in keepers_df.groupby("SampleID"):
    n_dup = sdf.shape[0]
    if n_dup > 1:
        keep = sdf.sort_values(by =['GISAID_Accession_Number','coverage_breadth','sequencing_depth_avg'], 
                               ascending=[False,False,False], na_position='last').iloc[0]
        print("  {:<8}  {:<8}  {:<10}  {:<4.1f}".format(n, n_dup, keep["SampleID"], keep["sequencing_depth_avg"]))
    else:
        keep = sdf.iloc[0]
    l_dfs.append(keep)
keepers_df = pd.concat(l_dfs, 1).transpose()
print("  Submittable samples: %d" % keepers_df.shape[0])
print("Done.")
print("")

print("Highlighting submittable samples")
keeperslist_df = keepers_df.filter(['SeqID','SeqRun','SeqBarcode'])
keeperslist_df["Submittable"] = True
alldata_df = pd.merge(left=merged_df,
                     right=keeperslist_df,
                     left_on=["SeqRun", "SeqBarcode", "SeqID"],
                     right_on=["SeqRun", "SeqBarcode", "SeqID"],
                    how='outer')
alldata_df = alldata_df.fillna({'Submittable' : False})
alldata_df = alldata_df.drop('ExcludeSample', axis=1)

#WRITE RESULTS
print("Writing all sequencing data...")
output_fn = "allsequencedata.tsv"
alldata_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = '\t', index=False)
print("  To: %s" % os.path.join(gisaid_dir, output_fn))
print("Done.")
print("")

#Calculate per run summaries
print("-" * 80)
print("Summarising output per run...")
print("  For all samples added to a sequencing run")
sample_sum = seqsamples_df.groupby("SeqRun").agg(
    samples=('SeqID', 'count')
)

print("  For all samples with a consensus")
sequenced_sum = sequenced_df.groupby('SeqRun').agg(
    all_consensus=('SeqID','count'),
    all_depth_min=('assembly_coverage_depth',min),
    all_depth_max=('assembly_coverage_depth',max),
    all_score_mean=('qc.overallScore','mean')
)

print("  For all samples with depth >{:<2}X and breadth >{:<2}% ".format(qc_depth, qc_breadth))
qc = sequenced_df[(sequenced_df.sequencing_depth_avg > qc_depth) & 
                    (sequenced_df.coverage_breadth > qc_breadth)]
qc_sum = qc.groupby('SeqRun').agg(
    qc_consensus=('SeqID','count'),
    qc_depth_min=('assembly_coverage_depth',min),
    qc_depth_max=('assembly_coverage_depth',max),
    qc_score_mean=('qc.overallScore','mean')
)

#Merge datasets together
sum1_df = pd.merge(left=sample_sum,
                     right=sequenced_sum,
                     left_on=["SeqRun"],
                     right_on=["SeqRun"],
                     how="outer")
sum1_df['all_success'] = sum1_df['all_consensus']/sum1_df['samples']

runsummary_df = pd.merge(left=sum1_df,
                     right=qc_sum,
                     left_on=["SeqRun"],
                     right_on=["SeqRun"],
                     how="outer")
runsummary_df['qc_success'] = runsummary_df['qc_consensus']/runsummary_df['samples']

print("  Done")

# WRITE RESULTS
print("  Writing results...")
runsummary_fn = "runsummary.csv"
runsummary_df.to_csv(os.path.join(gisaid_dir, runsummary_fn), index=True)
print("  To: %s" % os.path.join(gisaid_dir, runsummary_fn))
print("Done.")
print("")

print("-" * 80)
print("Runtime: %s" % str(datetime.timedelta(seconds=time.time() - start_time)))
print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("=" * 80)

