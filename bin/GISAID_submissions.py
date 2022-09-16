#!/usr/bin/env python3
# Complete submissions for GISAID submission
# -------------------------------------------------------------

import os
import sys
import time
import datetime
import getopt
import pandas as pd
import numpy as np
from Bio import SeqIO

print("=" * 80)
print("Prepare GISAID submission files")
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
    opts, args = getopt.getopt(sys.argv[1:], ":d:s:")
except getopt.GetoptError:
    print("  Error parsing options.")
    sys.exit(2)

for opt, value in opts:
    if opt == "-d":
        seqdata_dir = value
    elif opt == "-s":
        submission_dir = value
    else:
        print("  Parameter %s not recognized." % opt)
        sys.exit(2)

#Ensure all required variables have been parsed    
try:
    seqdata_dir
except:              # Whatever your error is (mine was a keyError)
    print()
    print("-d flag missing for sequence data directory" )
    print("using default ~/Dropbox/Zambia_nCoV/Data/Sequencing_Data")
    seqdata_dir = os.path.expanduser('~/Dropbox/Zambia_nCoV/Data/Sequencing_Data')
    if os.path.isdir(seqdata_dir):
        print("  Sequence data directory found: %s" % seqdata_dir)
    else:
        print("  Unable to recognise sequence data directory: %s" % seqdata_dir)
        sys.exit(2)
try:
    submission_dir
except:              # Whatever your error is (mine was a keyError)
    print()
    print("-s flag missing for submission directory" )
    print("  using default ~/Dropbox/Zambia_nCoV/Data/Submissions")
    submission_dir = os.path.expanduser('~/Dropbox/Zambia_nCoV/Data/Submissions')
    if os.path.isdir(submission_dir):
        print("  Submissions directory found: %s" % submission_dir)
    else:
        print("  Unable to recognise submissions directory: %s" % submission_dir)
        sys.exit(2)
print()
print("Done.")

def FormatDates (df,datefield):
    Date = df[datefield][df[datefield].notnull()].iloc[0]
    print("   Date field example: %s" % Date)
    if Date[2] == '/':
        print("Date format appears to be DD/MM/YYYY (Excel)")
        df[datefield] = pd.to_datetime(df[datefield], format='%d/%m/%Y')
    elif Date[4] == '-':
        print("Date format appears to be YYYY-MM-DD (native output)")
        df[datefield] = pd.to_datetime(df[datefield], format='%Y-%m-%d')
    else:
        print("Unknown date format")

print("-" * 80)
samples_seq_fn = "allsequencedata.csv"
print("Reading in %s file" % samples_seq_fn)
samplemeta_df = pd.read_csv(os.path.join(seqdata_dir,samples_seq_fn))

#Format Date fields
FormatDates(df=samplemeta_df, datefield="SpecimenDate")
FormatDates(df=samplemeta_df, datefield="SeqDate")

print("-" * 80)
print("Filtering for samples to submit")
samplemeta_df = samplemeta_df[samplemeta_df['Status'] == 'Awaiting Submission']
samplemeta_df.reset_index(inplace=True, drop=True)

if samplemeta_df.shape[0] == 0 :
    print("-" * 80)
    print("Exiting script as no entries present for submission")
    print("-" * 80)
    print("Runtime: %s" % str(datetime.timedelta(seconds=time.time() - start_time)))
    print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 80)
    sys.exit(0)
print("%d samples identified" % samplemeta_df.shape[0])

print("-" * 80)
print("Generating data for submission")


#Pull out the year
samplemeta_df["Year"] = samplemeta_df['SpecimenDate'].dt.year.astype('Int64')

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
#gisaid_df["covv_patient_status"] = samplemeta_df['PatientStatus']
gisaid_df["covv_specimen"] = "Nasopharyngeal swab"
gisaid_df["covv_outbreak"] = ""
gisaid_df["covv_last_vaccinated"] = ""
gisaid_df["covv_treatment"] = ""
gisaid_df["covv_seq_technology"] = "Nanopore MinION"
gisaid_df["covv_assembly_method"] = "ARTIC Field Workflow"
gisaid_df["covv_coverage"] = samplemeta_df['sequencing_depth_avg'].astype('int')
gisaid_df["covv_orig_lab"] = "University of Zambia, School of Veterinary Medicine"
#gisaid_df["covv_orig_lab"] = "Right to Care Zambia"
gisaid_df["covv_orig_lab_addr"] = "University of Zambia, School of Veterinary Medicine, Gt East Road Campus, Lusaka, Zambia"
#gisaid_df["covv_orig_lab_addr"] = "Right to Care, Lusaka, Zambia"
gisaid_df["covv_provider_sample_id"] = samplemeta_df["SeqID"]
#gisaid_df["covv_provider_sample_id"] = samplemeta_df["Epid Number"]
gisaid_df["covv_subm_lab"] = "UNZAVET and PATH"
gisaid_df["covv_subm_lab_addr"] = "University of Zambia, School of Veterinary Medicine, Gt East Road Campus, Lusaka, Zambia"
gisaid_df["covv_subm_sample_id"] = samplemeta_df["SampleID"]
gisaid_df["covv_authors"] = "Mulenga Mwenda, Ngonda Saasa, Daniel Bridges, ZNPHI and ZGSC"
#gisaid_df["covv_authors"] = "Mulenga Mwenda, Ngonda Saasa, Mupila Z, Kwenda G, Gill CJ, Mwananyanda L, Daniel Bridges"
gisaid_df["covv_comment"] = ""
gisaid_df["comment_type"] =""
print("   Done")

#WRITE RESULTS
print("-" * 80)
print("  Writing out metadata submission file for all submittable sequences...")
output_fn = "GISAID_Submission_Data.csv"
gisaid_df.to_csv(os.path.join(submission_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(submission_dir, output_fn))
print("   Done.")
print("")

#Combine dfs and drop all unnecessary columns
print("-" * 80)
samplemeta_df = samplemeta_df[['SampleID','SeqID','Notes','lineage']]
gisaid_df = gisaid_df[['covv_virus_name','covv_subm_sample_id']]
translate_df = pd.merge(samplemeta_df, gisaid_df, how='inner', left_on='SampleID', right_on='covv_subm_sample_id')
translate_df = translate_df.drop(columns='covv_subm_sample_id')

print("  Writing out translation file for filtering fasta and replacing SeqID with virus name")
output_fn = "GISAID_Translate.csv"
translate_df.to_csv(os.path.join(submission_dir, output_fn), sep = ',', index=False)
print("  To: %s" % os.path.join(submission_dir, output_fn))
print("Done.")
print("")

translate_df = translate_df.drop(columns=['Notes','lineage'])
print("-" * 80)
print("  Generating filtered fasta file for GISAID submission")
#Generate a dictionary of translationn terms
translate_dict = translate_df.set_index('SeqID').to_dict(orient='index')
#Set-up variables
count = 0
include = 0
retain = []
fasta_fn = "allsequences.fasta"

for record in SeqIO.parse(os.path.join(seqdata_dir, fasta_fn), "fasta"):
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

print("-" * 80)
print("  Writing out fasta file for submission")
#Output files
subseq_fn = "GISAID_Submission.fasta"
SeqIO.write(retain, os.path.join(submission_dir,subseq_fn), "fasta")
print("  To: %s" % os.path.join(submission_dir, subseq_fn))
print("Done.")
print("") 

print("-" * 80)
print("Runtime: %s" % str(datetime.timedelta(seconds=time.time() - start_time)))
print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("=" * 80)
