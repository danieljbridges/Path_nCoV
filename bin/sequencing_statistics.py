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
import time
import datetime
import getopt
import pandas as pd
import numpy as np
from Bio import SeqIO
from gisaid import *

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

# MERGE WITH SAMPLE LIST
print("-" * 80)
print("Merging sequenced samples with sample list...")
print("  No. samples...")
print("    ...in sample list: %d" % seqsamples_df.shape[0])
print("    ...with consensus sequence: %d" % stats_df.shape[0])
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

print("-" * 80)
print("Isolating submittable samples from controls, unsequenced samples, those marked for exclusion or failing qc parameters")
keepers_df = merged_df.query("sequencing_depth_avg >= @qc_depth" +
                             "& coverage_breadth >= @qc_breadth" +
                             "& ExcludeSample != 'Y'" +
                             "& Type == 'Sample'" ,
                             inplace=False)
print("    Samples remaining: %d" % keepers_df.shape[0])
print("Done")
print("")

#Highlight duplicates
print("Identifying highest depth concensus where samples have been sequenced multiple times:")

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
print("")
print("  Total Submittable samples: %d" % keepers_df.shape[0])
print("Done.")
print("")

print("Marking submittable samples in df")
keeperslist_df = keepers_df.filter(['SeqID','SeqRun','SeqBarcode'])
keeperslist_df["Submittable"] = True
alldata_df = pd.merge(left=merged_df,
                     right=keeperslist_df,
                     left_on=["SeqRun", "SeqBarcode", "SeqID"],
                     right_on=["SeqRun", "SeqBarcode", "SeqID"],
                    how='outer')
alldata_df = alldata_df.fillna({'Submittable' : False})

#WRITE RESULTS
print("Writing all sequencing data...")
output_fn = "allsequencedata.csv"
alldata_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(gisaid_dir, output_fn))
print("Done.")
print("")

#Merging in metadata
print("-" * 80)
print("Incorporating metadata...")

metadata_fn="Metadata.csv"
metadata_df = pd.read_csv(os.path.join(rampart_dir, metadata_fn))

#Add in column highlighting any missing data
check_cols = ["Province", "District", "SpecimenDate"] #List of key columns
missing_data = [] #To house the new column of data
for _, row in metadata_df.iterrows(): #Iterate over rows
    m = [] #List to hold output from an individual column
    for col in check_cols: #Iterate over columns we want to check
        if row[col] != row[col]: #Only NaN is not equal to itself
            m.append(col) #add data to list
    missing_data.append(", ".join(m)) #Join as a new string
metadata_df["MissingMetadata"] = missing_data

print("  Metadata identified for : %d samples" % metadata_df.shape[0])
print("  Sequencing data for : %d samples" % keepers_df.shape[0])
print("  Merging metadata with sequence data...")
samplemeta_df = pd.merge(left=metadata_df,
                     right=keepers_df,
                     left_on=["SampleID"],
                     right_on=["SampleID"],
                    how='inner')
print("Done")
print("  Total samples retained: %d" % keepers_df.shape[0])

#Change date fields from str to date
samplemeta_df['SeqDate'] = pd.to_datetime(samplemeta_df['SeqDate'], format='%d/%m/%Y')
samplemeta_df['SpecimenDate'] = pd.to_datetime(samplemeta_df['SpecimenDate'], format='%d/%m/%Y')
#Sanity check on sample date and sequencing date
samplemeta_df['DateCheck'] = samplemeta_df['SeqDate'] < samplemeta_df['SpecimenDate']
if samplemeta_df['DateCheck'].unique().shape[0] > 1:
    print("ERROR: %d records have a SeqDate before the SpecimenDate" % samplemeta_df[samplemeta_df.DateCheck==False].shape[0])
else:
    print("All records have a SeqDate after the SpecimenDate")
    samplemeta_df.drop(columns = ['DateCheck'], inplace = True)
print("")

#WRITE RESULTS
print("  Writing out all sequence data with metadata...")
output_fn = "Samples_Sequenced_With_Metadata.csv"
samplemeta_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(gisaid_dir, output_fn))
print("Done.")
print("")

print("-" * 80)
print("Generating GISAID submission files")

print("%d samples identified" % samplemeta_df.shape[0])
#Drop all entries without a date
samplemeta_df.query("SpecimenDate > datetime.datetime(2000,1,1)",inplace = True)
print("   %d samples with a date" % samplemeta_df.shape[0])
#Drop all entries previously submittedprint(l)
samplemeta_df.query("GISAID_Accession_Number != GISAID_Accession_Number",inplace = True)
print("   %d samples without an accession number" % samplemeta_df.shape[0])
#Drop all entries without a location
samplemeta_df.dropna(axis=0, subset=['Province'], inplace=True)
print("   %d samples with a location" % samplemeta_df.shape[0])
#Pull out the year
samplemeta_df["Year"] = samplemeta_df['SpecimenDate'].dt.year.astype('Int64')

#Reset the index
samplemeta_df.reset_index(inplace=True, drop=True)
#Generate a Province / District location
l = []
for _, row in samplemeta_df[['Province','District']].iterrows():
    l.append(pd.Series(row).str.cat(sep='/'))
samplemeta_df["Location"] = pd.DataFrame(l)

#Fill in NA entries and replace keys with correct values
samplemeta_df['Sex'].fillna(value="Unknown", inplace=True)
samplemeta_df['Sex'].replace({"M":"Male","F":"Female"}, inplace = True)
samplemeta_df['Age'].fillna(value="Unknown", inplace=True) 
samplemeta_df['PatientStatus'].fillna(value="Unknown", inplace=True)
#Reset the index so you can add rows correctly into the new dataframe
samplemeta_df.reset_index(inplace=True, drop=True)

#Create an empty dataframe with length of samplemeta_df
gisaid_df = pd.DataFrame(index=np.arange(samplemeta_df.shape[0]), columns=np.arange(0))
# pd.DataFrame(index=np.arange(1), columns=np.arange(8))
gisaid_df["submitter"] = "djbridges"
gisaid_df["fn"] = "GISAID_Submission_Data.csv"
gisaid_df["covv_virus_name"] = "hCoV-19/Zambia/ZMB-"+ samplemeta_df['SampleID'].astype('str') + "/" + samplemeta_df['Year'].astype('str')
gisaid_df["covv_type"] = "betacoronavirus"
gisaid_df["covv_passage"] = "Original"
gisaid_df["covv_collection_date"] = samplemeta_df['SpecimenDate']
gisaid_df["covv_location"] = "Africa/Zambia/" + samplemeta_df['Location']
gisaid_df["covv_add_location"] = ""
gisaid_df["covv_host"] = "Human"
gisaid_df["covv_add_host_info"] = ""
gisaid_df["covv_gender"] = samplemeta_df['Sex']
gisaid_df["covv_patient_age"] = samplemeta_df['Age']
gisaid_df["covv_patient_status"] = "Unknown"
gisaid_df["covv_specimen"] = "Nasopharyngeal swab"
gisaid_df["covv_outbreak"] = ""
gisaid_df["covv_last_vaccinated"] = ""
gisaid_df["covv_treatment"] = ""
gisaid_df["covv_seq_technology"] = "Nanopore MinION"
gisaid_df["covv_assembly_method"] = "ARTIC Field Workflow"
gisaid_df["covv_coverage"] = samplemeta_df['sequencing_depth_avg'].astype('int')
gisaid_df["covv_orig_lab"] = "University of Zambia, School of Veterinary Medicine"
gisaid_df["covv_orig_lab_addr"] = "University of Zambia, School of Veterinary Medicine, Gt East Road Campus, Lusaka, Zambia"
gisaid_df["covv_provider_sample_id"] = samplemeta_df["SeqID"]
gisaid_df["covv_subm_lab"] = "UNZAVET and PATH"
gisaid_df["covv_subm_lab_addr"] = "University of Zambia, School of Veterinary Medicine, Gt East Road Campus, Lusaka, Zambia"
gisaid_df["covv_subm_sample_id"] = samplemeta_df["SampleID"]
gisaid_df["covv_authors"] = "Mulenga Mwenda-Chimfwembe, Ngonda Saasa, Daniel Bridges, ZNPHI and ZGSC"
gisaid_df["covv_comment"] = ""
gisaid_df["comment_type"] =""
print("Done")

#WRITE RESULTS
print("  Writing out gisaid metadata submission file for all submittable sequences...")
output_fn = "GISAID_Submission_Data.csv"
gisaid_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(gisaid_dir, output_fn))
print("Done.")
print("")

#Combine dfs and drop all unnecessary columns
samplemeta_df = samplemeta_df[['SampleID','SeqID']]
gisaid_df = gisaid_df[['covv_virus_name','covv_subm_sample_id']]
translate_df = pd.merge(samplemeta_df, gisaid_df, how='inner', left_on='SampleID', right_on='covv_subm_sample_id')
translate_df = translate_df.drop(columns='covv_subm_sample_id')

print("  Writing out translation file for filtering fasta and replacing SeqID with virus name")
output_fn = "GISAID_Translate.csv"
translate_df.to_csv(os.path.join(gisaid_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(gisaid_dir, output_fn))
print("Done.")
print("")

print("  Generating filtered fasta file for GISAID submission")
#Generate a dictionary of translationn terms
translate_dict = translate_df.set_index('SeqID').to_dict(orient='index')
#Set-up variables
count = 0
include = 0
retain = []
fasta_fn = "allsequences.fasta"

for record in SeqIO.parse(os.path.join(gisaid_dir, fasta_fn), "fasta"):
    count += 1
    if translate_dict.get(record.name): #If record.name is in the dictionary
        #Pull out virus name given
        vname = translate_dict.get(record.name).get('covv_virus_name')
        #Update the record id and description otherwise if not the same they are concat in fasta header
        record.id = vname
        record.description = vname
        #Add to list for export
        retain.append(record)
        include += 1
        
print("  %d sequences identified" % count)
print("    %d sequences retained" % include)
#Output files
allseqs_fn = "GISAID_Submission.fasta"
SeqIO.write(retain, os.path.join(data_dir,allseqs_fn), "fasta")
print("All %d sequences output to %s" % (len(retain),allseqs_fn))
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

