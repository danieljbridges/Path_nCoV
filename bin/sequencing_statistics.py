#!/usr/bin/env python
# coding: utf-8
# ----------------------------------------
# Compute statistics for GISAID submission
# ----------------------------------------
#
# Example usage:
#   python run_gisaid-statistics.py -d data/WGS
#
# JHendry and DBridges


import os
import sys
import time
import datetime
import getopt
import pandas as pd
import numpy as np
from Bio import SeqIO
from gisaid import *
import configparser

print("=" * 80)
print("PATH-UNZAVET bioinformatic pipeline to calculate / combine sequencing statistics")
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

#Loading configuration file
print("=" * 80)
config_file = os.path.join("../configs/", "default.ini")
config = configparser.ConfigParser()
config['DEFAULT'] = {'qc_breadth': 50,
                     'qc_depth': 50}
if not os.path.isfile(config_file):
    print ("Missing configuration file - using defaults")
    qc_depth = int(config['DEFAULT']['qc_depth'])
    qc_breadth = int(config['DEFAULT']['qc_breadth'])
else:
    print ("Found configuration file")
    config.read(config_file)
    qc_depth = int(config['SubmissionCriteria']['qc_depth'])
    qc_breadth = int(config['SubmissionCriteria']['qc_breadth'])
print ("   Using qc_depth cutoff of %d" % qc_depth)
print ("   Using qc_breadth cutoff of %d" % qc_breadth)
print("")

# EXAMINE CONTENTS
print("=" * 80)
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
dts = [] #data per sample
dtr = [] #data per run

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
        
        #Dictionary for a count of the number of reads that are unclassified
        stats_run = {}
        stats_run["SeqRun"] = r
        
        #Need to count the number of unclassified reads just once per run
        unclass_dir = os.path.join(d.replace("processed", "fastq/unclassified"))
        t = 0 #count number of reads
        try:
            for h, u in enumerate(os.listdir(unclass_dir)): #for each fastq
                current = calc_fastq_total_reads(os.path.join(unclass_dir, u))
                t = t+current
            stats_run.update({"total_reads": t})
        except:
            stats_run.update({"total_reads": 0})
        
        # Store unclass results
        dtr.append(stats_run)
        
        #Iterate over samples
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
            
            # Identify fastq files for analysis
            b_fn = os.path.join(d.replace("processed", "fastq"), "C%02d_barcode%02d.fastq" % (int(r[1:]), b))
            try:
                # Calc. sequencing depth from .fastq
                sequencing_depth_avg_fastq = calc_avg_seq_depth(b_fn, genome_length=stats_dt["ref_genome_length"])
                stats_dt.update({"sequencing_depth_avg_fastq": sequencing_depth_avg_fastq})
                # Calc. number of reads per barcode
                total_reads = calc_fastq_total_reads(b_fn)
                stats_dt.update({"total_reads": total_reads})
            except:
                stats_dt.update({"sequencing_depth_avg_fastq": 0})
                stats_dt.update({"total_reads": 0})
            
            #Calculate coverage from fasta file
            
            # Store
            dts.append(stats_dt)
        
        sys.stdout.flush()
        sys.stdout.write("\n")
        
# Create data frame
stats_df = pd.DataFrame(dts)
print("Done.")

fa_path = os.path.join(gisaid_dir, "allsequences.fasta")
if os.path.isfile(fa_path):    
    print("Calculating statistics from the multi-fasta file for all sequences")
    #Generate statistics from the multi-fasta file
    fasta_dt = fasta_stats(fa_path)
    #Turn dictionary into a dataframe in correct orientation
    fasta_df = pd.DataFrame.from_dict(fasta_dt, orient ='index')
    #Reset the index so it is a column you can match on
    fasta_df.reset_index(inplace=True)
    print("   Stats generated for %d sequences" % len(fasta_df))

    #Merge with the other gisaid stats
    print("Merging fasta stats with bam / fastq stats (stats_df)")
    print("   %d records in stats_df" % len(stats_df))
    stats_df = pd.merge(stats_df, fasta_df, how='outer', left_on='SeqID', right_on='index')
    print("   %d records remaining" % len(stats_df))
    stats_df.drop(columns = ['index'], inplace = True)


#Check that all the sequence IDs are unique
if len(stats_df["SeqID"]) == len(set(stats_df["SeqID"])):
    print("    All sequence IDs are unique")
else:
    print("    Duplicate sequence IDs identified")
    l_func = lambda x, y: list((set(x)- set(y))) + list((set(y)- set(x))) 
    non_match = l_func(stats_df["SeqID"], stats_df["SeqID"].unique()) 
    print("    Non-match elements: ", non_match)
    print("Exiting script")
    exit()
print("")
print("Statistics for a total of %s samples have been calculated" % stats_df.shape[0])
stats_fn = "intermediates/stats.csv"
stats_df.to_csv(os.path.join(gisaid_dir, stats_fn), index=False)
print("    Data output to %s" % stats_fn)
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
non_match_gp = l_func(stats_df["SeqID"], pango_df["taxon"]) 
non_match_gn = l_func(stats_df["SeqID"], nextclade_df["seqName"]) 
if len(non_match_gp) > 0 :
    print("  A total of %s sequence ids do not match between gisaid and PANGO" % len(non_match_gp))
    #exit()
elif len(non_match_gn) > 0 :
    print("  A total of %s sequence ids do not match between gisaid and nextclade" % len(non_match_gn))
    #exit()
else:
    print("  All sequences match across the three files")

print("  Merging gisaid (n=%d), pango (n=%d) and nextclade (n=%d) data together" % (stats_df.shape[0], pango_df.shape[0], nextclade_df.shape[0]))
Merge1 = pd.merge(stats_df, pango_df, how='outer', left_on='SeqID', right_on='taxon')
Merge2 = pd.merge(Merge1, nextclade_df, how='outer', left_on='SeqID', right_on='seqName')
print("  Total remaining records = %d" % (Merge2.shape[0]))
print("  Removing unwanted column headings from final summary file of samples sequenced")
gisaidcols = ['ref_genome_length']
pangocols = (['taxon'])
nextcladecols = ['seqName']
dropcols = gisaidcols + pangocols + nextcladecols
sequenced_df = Merge2.drop(dropcols, axis=1)
print("Done.")
print("")

#Read in and edit metadata
print("-" * 80)
print("Reading in metadata...")

metadata_fn="Metadata.csv"
metadata_df = pd.read_csv(os.path.join(rampart_dir, metadata_fn))

# MERGE ALL DATA SOURCES
print("-" * 80)
print("Merging Sample data, SequenceData, and Metadata")
print("-" * 80)
print("   %d sample records" % samples_df.shape[0])
print("   %d sequence data records" % sequenced_df.shape[0])
print("   %d metadata records" % metadata_df.shape[0])

# Retain all records
alldata_df = pd.merge(left=samples_df,
                     right=sequenced_df,
                     left_on=["SeqRun", "SeqBarcode", "SeqID"],
                     right_on=["SeqRun", "SeqBarcode", "SeqID"],
                    how='outer')
# Only retain records that match with what was sequenced
alldata_df = pd.merge(left=alldata_df,
                     right=metadata_df,
                     left_on=["SampleID"],
                     right_on=["SampleID"],
                    how='left')

print("%d all_data records" % alldata_df.shape[0])
if samples_df.shape[0] != alldata_df.shape[0]:
    print("ERROR: there should be %d records" % samples_df.shape[0])
print("-" * 80)

#Highlight records with missing metadata
print("Highlight missing metadata fields")
check_cols = ["Province", "District", "SpecimenDate"] #List of key columns
missing_data = [] #To house the new column of data
for _, row in alldata_df.iterrows(): #Iterate over rows
    m = [] #List to hold output from an individual column
    for col in check_cols: #Iterate over columns we want to check
        if row[col] != row[col]: #Only NaN is not equal to itself
            m.append(col) #add data to list
    missing_data.append(", ".join(m)) #Join as a new string
alldata_df["MissingMetadata"] = missing_data
print ("Done")
print("-" * 80)

#Highlight duplicates
print("Identifying highest depth consensus where samples have been sequenced multiple times:")
print("")

l_dfs = []
print("  {:<8}  {:<8}  {:<15}  {:<4}  {:<4}  {:<10}  {:<4}  {:<4}  {:<4}"
      .format("Sample", "No. dup.", "SeqID", "Breadth", "Depth", "Note", "Top_SeqID","Top_Breadth", "Top_Depth",))

for n, sdf in alldata_df.groupby("SampleID"):
    n_dup = sdf.shape[0]
    if n_dup > 1:
        keep = sdf.sort_values(by =['GISAID_Accession_Number','coverage_breadth','sequencing_depth_avg'], 
                               ascending=[False,False,False], na_position='last').iloc[0]
        #Identify best SeqID scoring entry
        top = sdf.sort_values(by =['coverage_breadth','sequencing_depth_avg'], 
                                   ascending=[False,False], na_position='last').iloc[0]
        #Check it correlates to the one selected to keep and warn if not else print summary
        if top['SeqID'] != keep['SeqID'] :
            print("  {:<8}  {:<8}  {:<15}  {:<4.0f}     {:<4.0f}     {:<4}  {:<15}  {:<4.0f}     {:<4.0f}".format(
                n, n_dup, keep["SeqID"], keep["coverage_breadth"], keep["sequencing_depth_avg"], 
                "WARNING", top['SeqID'], top['coverage_breadth'], top['sequencing_depth_avg']))
        else:
            print("  {:<8}  {:<8}  {:<15}  {:<4.0f}     {:<4.0f}".format(
                    n, n_dup, keep["SeqID"], keep["coverage_breadth"], keep["sequencing_depth_avg"]))
    else:
        keep = sdf.iloc[0]
    l_dfs.append(keep)

topscores_df = pd.concat(l_dfs, 1).transpose()

print("-" * 80)
print("Marking top scoring records for each SampleID")
alldata_df["TopScore"] = np.where(alldata_df["RecordID"].isin(topscores_df["RecordID"]),True,False)

print("Marking records where sequence sample passes QC standards (breadth %d%%, depth %dx)" % (qc_breadth, qc_depth))
alldata_df["QC_Sequence"] = np.where((alldata_df["sequencing_depth_avg"] > qc_depth) &
                                     (alldata_df["coverage_breadth"] >= qc_breadth),True,False)

print("Marking status of each row")
alldata_df.loc[alldata_df['GISAID_Accession_Number'].notnull(), 'Status'] = 'Published'
alldata_df.loc[alldata_df['Type'] == 'Control', 'Status'] = 'Controls'
alldata_df.loc[alldata_df['ExcludeSample']=='Y', 'Status'] = 'Excluded'
alldata_df.loc[(alldata_df['Status'].isnull()) &
               (alldata_df['QC_Sequence']==True) &
               (alldata_df['TopScore'] == True) &
               (alldata_df['SpecimenDate'].notnull()), 'Status'] = 'Awaiting Submission'
alldata_df.loc[(alldata_df['Status'].isnull()) & 
               (alldata_df['TopScore'] == True) &
               (alldata_df['QC_Sequence']==True) &
               (alldata_df['SpecimenDate'].isnull()), 'Status'] = 'Requires metadata for submission'
alldata_df.loc[(alldata_df['Status'].isnull()) &
               (alldata_df['SeqRun'].notnull()) , 'Status'] = 'Sequenced'
alldata_df.loc[alldata_df['Status'].isnull() , 'Status'] = 'PCR amplified'

print("-" * 80)

print("Tidying up various data fields")
#Change date fields from str to date
alldata_df['SeqDate'] = pd.to_datetime(alldata_df['SeqDate'], format='%d/%m/%Y')

#Drop date field as excel keeps changing format and rebuild from separate cols
alldata_df.drop(labels='SpecimenDate', axis =1)
alldata_df['SpecimenDate'] = pd.to_datetime(alldata_df['SpecimenDate'], format='%d/%m/%Y')
#Sanity check on sample date and sequencing date
alldata_df['DateError'] = alldata_df['SeqDate'] < alldata_df['SpecimenDate']

if alldata_df[alldata_df.DateError==True].shape[0] > 0 :
    print("ERROR: %d records have a SeqDate before the SpecimenDate" % alldata_df[alldata_df.DateError==True].shape[0])
    print(alldata_df[alldata_df.DateError==True][['SeqID','SeqRun','SpecimenDate','SeqDate']])
else:
    print("All records have a SeqDate after the SpecimenDate")
    alldata_df.drop(columns = ['DateError'], inplace = True)
print("")

#WRITE RESULTS
print("Writing out all sequencing data...")
output_fn = "allsequencedata.csv"
alldata_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = ',', index=False)
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

#Create dataframes for barcoded and unclassified reads per run
bc = stats_df[['SeqRun','total_reads']].groupby('SeqRun').sum()
bc.rename(columns = {'total_reads':'barcoded_reads'}, inplace = True)
un = pd.DataFrame(dtr).groupby('SeqRun').sum()
un.rename(columns = {'total_reads':'unclassified_reads'}, inplace = True)
#Merge into runsummary_df
runsummary_df = pd.merge(left=runsummary_df,
                     right=bc,
                     left_on=["SeqRun"],
                     right_on=["SeqRun"],
                     how="outer")
runsummary_df = pd.merge(left=runsummary_df,
                     right=un,
                     left_on=["SeqRun"],
                     right_on=["SeqRun"],
                     how="outer")
#Calculate barcoding efficiency
runsummary_df["barcoding_efficiency"] = runsummary_df['barcoded_reads'] / (runsummary_df['unclassified_reads'] + runsummary_df['barcoded_reads'])

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

