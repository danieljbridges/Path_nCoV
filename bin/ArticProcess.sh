#!/bin/bash
# ----------------------------------------------------------------------
# Script to automate the artic processing pipeline for SARS-CoV-2 WGS
# Written by Daniel Bridges (danieljbridges@gmail.com)
#
#Changelog
# ----------------------------------------------------------------------
set -e #exit whenever a command exits with a non zero status
set -u #treat undefined variables as errors
set -o pipefail #pipe will be considered successful if all the commands are executed without errors

VERSION="0.3.4"
VERDATE="2021/02/25"
#ANSI escape codes: 
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
LG='\033[1;32m'
RED='\033[0;31m'
GREEN='\033[0;32m' 
BLUE='\033[0;34m' 
ORANGE='\033[0;33m' 
NC='\033[0m' # No Color

#==============FUNCTIONS==================

function Help {
   # Display Help
   printf "${RED}##### HELP NOTES: #####${NC}
   This script chains together the various processes of the artic network bioinformatic pipeline for each sample. It assumes the following data folder structure
   
   basefolder/1_Raw                         Raw data output from minion
   basefolder/2_SampleList_and_Rampart      List of all samples
   basefolder/3_Artic_Output                Processed files from artic pipeline
   basefolder/4_Consensus                   Consensus sequences multi-fasta per run
   basefolder/5_GISAID                      Sequences for uploading to GISAID
   
   Syntax: $(basename "$0") [-123|b|n|p|r|]
   Version: $VERSION, $VERDATE
   
    Mandatory args:
    -b      Basefolder for all inputs and outputs - please use an absolute path
    -r      Runname to prefix output files, folder names etc

    Optional args:
    -h      Print this help
    -m      Process using the medaka pipeline rather than with nanopolish
    -p      Primers used during PCR [Sanger | Artic]
    -s      Location directory of the primer scheme used
    -v      Version

    Optional steps to omit:
    -1      Omit Step 1 (guppy_barcoder)
    -2      Omit Step 2 (guppyplex)
    -3      Omit Step 3 (artic minion)
    -4      Omit Step 4 (Consensus sequences)
    -5      Omit Step 5 (Sequencing statistics)
    -6      Omit Step 6 (Generate JSON file for RAMPART)
   "
}

function present {
    if [ $2 = 'd' ] ; then
        if [ ! -d $1 ] ; then
            printf "\n${RED}ERROR:${NC} $1 is not present\nExiting script.\n"
            exit
        fi
    elif [ $2 = 'f' ] ; then
        if [ ! -f $1 ] ; then
            printf "\n${RED}ERROR:${NC} $1 is not present \nExiting script.\n"
            exit
        fi
    fi
}

function check_mkdir {
    if [ ! -d $1 ] ; then
        echo -e "Making directory $1"
        mkdir -p $1
    fi
}

function check_var {
    if [ -z $1 ] || [ $1 = 1 ]; then
        Help
        printf "\n${RED}ERROR:${NC} $2 variable not supplied. See Help message above.\n Exiting script\n\n"
        exit   
    fi
}

function check_package {
    if [ $(which $1 | wc -l) = 0 ] ; then
        printf "${RED}ERROR:${NC} $1 not present. $2 \n Exiting script\n\n"
        exit
    fi
}

#==============Recognise any CLI options==================
while getopts 'vhp:b:r:s:m123456' opt; do
  case "$opt" in
    v)
      printf "${ORANGE}ArticProcess.sh Version $VERSION, $VERDATE ${NC}\n"
      exit
      ;;
    b)
      BASEFOLDER=$OPTARG
      ;;
    p)
      PRIMERS=$OPTARG
      ;;
    s)
      PRIMERDIR=$OPTARG
      ;;
    h)
      Help
      exit
      ;;
    m)
      MEDAKA=1
      ;;
    1)
      S1=0
      ;;
    2)
      S2=0
      ;;
    3)
      S3=0
      ;;
    4)
      S4=0
      ;;
    5)
      S5=0
      ;;
    6)
      S6=0
      ;;
    r)
      RUNNAME=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      #echo "script usage: $(basename $0) [-l] [-h] [-a somevalue]" >&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

#==============CHECK ALL VARIABLES AND REQUIREMENTS ARE MET==================
#Create a space between cmd
printf "\n${ORANGE}##### Running ArticProcess.sh Version $VERSION ##### ${NC}\n"

#################################################
#Enter default values if parameter not defined or build from exisiting args
# These are MANDATORY

check_var ${BASEFOLDER:- 1} "-b (basefolder)"
BASEFOLDER=`realpath $BASEFOLDER`
present $BASEFOLDER "d"

check_var ${RUNNAME:- 1} "-r (runname)"

#Fill in all undefined variables
S1=${S1:- 1}
S2=${S2:- 1}
S3=${S3:- 1}
S4=${S4:- 1}
S5=${S5:- 1}
S6=${S6:- 1}
MEDAKA=${MEDAKA:- 0}            # 3
PRIMERDIR=${PRIMERDIR:- 1}      # 3
PRIMERS=${PRIMERS:- 1}          # 2, 3
ARTIC_OUT="$BASEFOLDER/3_Artic_Output/$RUNNAME"
RAWDATADIR="$BASEFOLDER/1_Raw/$RUNNAME"


###############
#                   Needed by step
# Packages
###### guppy_barcoder    1
# artic             1,2,3

# files
# Sequencing summary  3  

# CLI VARIABLES
#      BASEFOLDER   1-6if (( ${#SAMPLES[@]} == 0 )); then
#      PRIMERS      2, 3
#      PRIMERDIR    3
#      RUNNAME      1-6

# FASTQRAW          1
# PRIMERS           2, 3
# MAX               2, 3
# MIN               2, 3
# PRIMERSCHEME      2, 3
# MEDAKA            3
# PRIMERDIR         3
# SAMPLEFILE        3
# SAMPLES           3
# BARCODES          3
# BASEFOLDER        1,2,3,4,5,6
# RUNNAME
#which:
# guppy_barcoder    1


#################################################
#Check necessary packages and folders for all functions
check_package "artic environment" "Please activate the appropriate environment e.g.:\n\n conda activate artic\n"

#################################################
#Determine if the necessary files and directories and arguments etc exist or create them if not

#Make directories for the demultiplexed, combined and size selected fast_q files

check_mkdir "$ARTIC_OUT/fastq"
check_mkdir "$ARTIC_OUT/processed"

present $RAWDATADIR "d"


###############
#Guppy Barcoder (step 1)
if [ $S1 = 1 ] ; then
    check_package "guppy_barcoder"
    FASTQRAW="$RAWDATADIR/fastq_pass"
    present $FASTQRAW "d"
fi

###############
#Guppy Plex (step 2) and step 3 identical
if [ $S2 = 1 ] || [ $S3 = 1 ]; then
    #Check primers and determine max and min sizes and primer scheme for analysis
    check_var $PRIMERS "-p (primers)"
    if [ $PRIMERS = "Sanger" ] ; then
        MAX=1250
        MIN=750
        PRIMERSCHEME=SARS-CoV-2/V1
    elif [ $PRIMERS = "Artic" ] ; then
        MAX=700
        MIN=400
        PRIMERSCHEME=SARS-CoV-2/V3
    else
        Help
        printf "${RED}ERROR:${NC} Unrecognised primer scheme\n\n"
        exit
    fi
fi

###############
#Artic processing (step 3)
if [ $S3 = 1 ] ; then
    SEQUENCINGSUMMARY=$(find $RAWDATADIR -name 'sequencing_summary*' -not -name '*.tmp')
    check_var ${SEQUENCINGSUMMARY:- 1} "Sequencing Summary"
    
    #Check primer scheme supplied
    check_var $PRIMERDIR "-s (primer scheme location)"
    PRIMERDIR=`realpath $PRIMERDIR`
fi

###############
#Artic processing and JSON  (step 3 and 6)
if [ $S3 = 1 ] || [ $S6 = 1 ] ; then
    #Identify csv file
    SAMPLEFILE="$BASEFOLDER/2_SampleList_and_Rampart/Samples_Sequenced.csv"
    present $SAMPLEFILE "f"
    #Check if all sample IDs are unique in the samples file
    readarray -t IDS < <(awk -F"," 'NR>1 {print$2}' $SAMPLEFILE | sort | uniq -d)
    printf '%s\n' "${IDS[@]}"
    if [ ${#IDS[@]} -gt 0 ]; then
        printf "${RED}ERROR:${NC} The following duplicate samples exist in $SAMPLEFILE. 
        Please correct so that every uniqueID (column 2) is unique to the entire list.\n"
        printf '%s\n' "${IDS[@]}"
        exit
    fi
    
    #Read in Samples
    readarray -t SAMPLES < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$2}' $SAMPLEFILE)
    if (( ${#SAMPLES[@]} == 0 )); then
        printf "\n${RED}ERROR:${NC} SAMPLES array is empty \nExiting script\n\n"
        exit
    fi
    
    #Read in barcodes
    readarray -t BARCODES < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$1}' $SAMPLEFILE)
    if (( ${#BARCODES[@]} == 0 )); then
        printf "\n${RED}ERROR:${NC} BARCODES array is empty \nExiting script\n\n"
        exit
    fi
    
    #Check for barcode duplicates within the run
    if [ `printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1' | wc -l` -gt 0 ] ; then
        printf "${RED}ERROR:${NC} The following duplicate barcodes are present in $SAMPLEFILE file for $RUNNAME run:\n"
        printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1'
        exit
    fi
fi
#Define logfile output
LOGFOLDER="$BASEFOLDER/Logs/$RUNNAME"
check_mkdir $LOGFOLDER
LOG=$(echo "$LOGFOLDER/${RUNNAME}_Step_")

printf "${GREEN}CHECKED:${NC}All required programs, files and locations are present.
${GREEN}CHECKED:${NC}No errors in sample list.\n"

#==============THE SCRIPT==================Se   
#STEP 1: Run the guppy barcoder to demultiplex into separate barcodes
if [ $S1 = 1 ] ; then
    printf "\n###### ${BLUE}Step 1: Running the guppy_barcoder to demultiplex the FASTQ files.${NC} ######\n\n"
    printf "${LG}guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files barcode_arrs_nb96.cfg${NC}\n" | tee "${LOG}1.log"
    
    guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files "barcode_arrs_nb96.cfg" | tee -a "${LOG}1.log"
    printf "\n###### ${GREEN}Step 1: guppy_barcoder completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 1: Skipping guppy_barcoder step${NC} ######\n\n"
fi

#STEP 2: Combine all identical barcode reads into same fasta and exclude sizes
if [ $S2 = 1 ] ; then
    printf "\n###### ${BLUE}Step 2: Combining demultiplexed files into a single fastq and excluding based on size.${NC} ######\n\n" | tee "${LOG}2.log"
    readarray -t S2DIRS < <(find $ARTIC_OUT/fastq -type d -name 'barcode[0-9]*')
    
    #Change directory as unable to redirect output from guppyplex
    cd $ARTIC_OUT/fastq
    
    #Run through the demuxed barcodes to combine into single fastq files
    for S2DIR in "${S2DIRS[@]}" ; do
        printf "${LG}artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME ${NC}\n" | tee -a "${LOG}2.log"
       artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME | tee -a "${LOG}2.log"
    done
    echo -e "\n###### ${GREEN}Step 2: guppyplex completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 2: Skipping guppyplex step${NC} ######\n\n"
fi

#STEP 3: Processing with artic minion pipeline
if [ $S3 = 1 ] ; then
    printf "\n###### ${BLUE}Step 3: Importing samplenames and processing with artic minion command.${NC} ######\n\n" | tee "${LOG}3.log"
    
    #Change directory
    cd $ARTIC_OUT/fastq/
    #Set the count
    COUNT=0
    
    #Cycle through each barcode and look for the matching file
    for i in "${BARCODES[@]}"; do 
        #Identify correct samplenames / barcodes and then run the processing
        BARCODE=`printf "%02.f" ${BARCODES[$COUNT]}` #pad left using float otherwise 0 preceeding = octal number
        SAMPLENAME=${SAMPLES[$COUNT]}
        FILE="${RUNNAME}_barcode${BARCODE}.fastq"
        
        #Increment count
        (( COUNT += 1 ))
        
        #Check that there is a file for this barcode (NTC may not have any reads)
        if [ ! -f $FILE ] ; then
            printf "${RED}WARNING:${NC} $FILE is not present - moving to next sample.\n" | tee -a "${LOG}3.log"
        else
            #Ensure fastq file is not empty
            FASTQLENGTH=`wc -l $FILE | sed s/$FILE//`  
            
            echo -e "\n\n \n${ORANGE} Processing barcode number $BARCODE, Sample $SAMPLENAME from file $FILE \n${NC}" | tee -a "${LOG}3.log"
            
            #Remove files that have not got enough data
            if [ $FASTQLENGTH -gt 1000 ] ; then
                #Run processing scheme
                if [ $MEDAKA = 1 ] ; then
                    printf "\n${GREEN}Using Medaka pipeline${NC}\n"
                    printf "${LG}artic minion --medaka --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE $PRIMERSCHEME $SAMPLENAME ${NC}\n" | tee -a "${LOG}3.log"
                    artic minion --medaka --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${LOG}3.log"
                    # 2>&1 redirects stderr to stdout
                else
                    printf "\n${GREEN}Using Nanopolish pipeline${NC}\n"
                    printf "${LG}artic minion --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME${NC}\n" | tee -a "${LOG}3.log"
                    artic minion --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${LOG}3.log"
                fi
                
                #Make samplename subdirectory if required
                if [ ! -d $SAMPLENAME ] ; then
                    mkdir $SAMPLENAME
                    echo "Making subdirectory $SAMPLENAME"
                fi
                
                #Move output from what was produced into its own directory
                find ./ -type f -name "$SAMPLENAME*" | xargs -I '{}'  mv {} "$SAMPLENAME"/
                #Move directory to another level for clarity
                mv $ARTIC_OUT/fastq/$SAMPLENAME $ARTIC_OUT/processed/$SAMPLENAME
            else
                printf "\n\n${RED}ERROR:${NC} Too few reads (n = $FASTQLENGTH) in File $FILE (Barcode $BARCODE, Sample $SAMPLENAME).\nAborting processing this file\n" | tee -a "${LOG}3.log"
            fi
        fi
    done
    echo -e "\n###### ${GREEN}Step 3: artic minion completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 3: Skipping artic minion step${NC} ######\n\n"
fi

#STEP 4: Concatenate consensus sequences
if [ $S4 = 1 ] ; then
    printf "\n###### ${BLUE}Step 4: Concatenating consensus sequences from artic pipeline.${NC} ######\n\n"
    cd $BASEFOLDER
    #Ensure folder exists
    CONSENSUS="$BASEFOLDER/4_Consensus"
    check_mkdir $CONSENSUS
    cd $CONSENSUS
    #Concatenate all of the fasta consensus data
    find /$ARTIC_OUT/processed -name '*.consensus.fasta' | xargs cat > $RUNNAME.consensus.fasta } | tee "${LOG}4.log"
    
    echo -e "\n###### ${GREEN}Step 4: consensus sequences compiled. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 4: Skipping consensus sequences ${NC} ######\n\n"
fi


#STEP 5: Generate GISAID stats
if [ $S5 = 1 ] ; then
    printf "\n###### ${BLUE}Step 5: Generating statistics on sequencing and identify sequences for upload. ${NC} ######\n\n"
    #Move to basefolder to ensure output is put int the correct location
    cd $BASEFOLDER
    #Run 
    run_gisaid-statistics.py -d $BASEFOLDER
    
    #Filter the list of all sequences to only retain the correct ones
    #seqkit grep -n -f FASTAHeadersFilter.csv all.fasta -o Filtered.fasta

    printf "\n###### ${GREEN}Step 5: Sequencing statistics compiled. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 5: Skipping generating sequencing statistics${NC} ######\n\n"
fi

if [ $S6 = 1 ] ; then
    printf "\n###### ${BLUE}Step 6: Generating JSON files for Rampart. ${NC} ######\n\n"

    COUNT=0
    JSONFILE="$BASEFOLDER/2_SampleList_and_Rampart/${RUNNAME}_run_configuration.json"
    
    #Start the JSON File
    printf "{
        \"title\": \"Run $RUNNAME\",
        \"basecalledPath\": \"fastq_pass\",
        \"samples\": [" > $JSONFILE
        
    for BARCODE in "${BARCODES[@]}"; do 
        SAMPLE=${SAMPLES[$COUNT]}
        #Ensure that all barcodes are 2 digit
        if [ ${#BARCODE} -lt 2 ]; then
            BARCODE="0$BARCODE"
        fi
        printf "    {
        \"name\": \"$SAMPLE ($BARCODE)\",
        \"description\": \"\",
        \"barcodes\": [ \"NB$BARCODE\" ]" >> $JSONFILE
        
        if [ ${#BARCODES[@]} == $COUNT ]; then
            echo -e "    }" >> $JSONFILE #If last entry then there shouldn't be a comma
        else
            echo -e "    }," >> $JSONFILE 
        fi
        #Increment the count
        #Increment count
        (( COUNT += 1 )) #let COUNT=COUNT+1
        printf "${BLUE}ADDED:${NC} Sample = $SAMPLE, Barcode = $BARCODE\n"
    done
    #Finish off the JSON File
    printf "  ]
    }" >> $JSONFILE
    printf "\n###### ${GREEN}Step 6: Rampart config saved to $JSONFILE. ${NC} ######\n\n"
fi

exit

######################## NOTES ##########
#TO DO
#Retain only the sequences that are final
#1. Filter to final multifasta
seqkit grep -n -f FASTAHeadersFilter.csv all.fasta > Filtered.fasta
#2. Replace headers
seqkit replace -p "(.+)" -r '{kv}' -i -K -k FASTAHeadersReplacements.csv Filtered.fasta > Submit.fasta

#perform pangolin lineages on what was created?
    #GISAID=`find ./ -name '*gisaid.csv' -printf "%f\n"`
    #readarray -t KEEPERS < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$2}' $SAMPLEFILE)
    #readarray -t SAMPLES < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$2}' $SAMPLEFILE)
    cd

#Create Samplejson list from the csv file

cat all.fasta | sed s'/\//_/g' | seqkit replace -p "(.+)" -r '{kv}' -K -k Replacements.tsv | seqkit grep -r -p ^SARS > Submission.fa

#Bug catching - list all defined variables
( set -o posix ; set ) | less



