#!/bin/bash
# ----------------------------------------------------------------------
# Script to automate the artic processing pipeline for SARS-CoV-2 WGS
# Version 0.3 03/02/2021
# Written by Daniel Bridges (danieljbridges@gmail.com)
#
#Changelog
# ----------------------------------------------------------------------
set -e #exit whenever a command exits with a non zero status
set -u #treat undefined variables as errors
set -o pipefail #pipe will be considered successful if all the commands are executed without errors

#==============FUNCTIONS==================

function Help {
   # Display Help
   echo -e "This script chains together the various processes of the artic network bioinformatic pipeline for each sample. It assumes the following data folder structure
   
   basefolder/1_Raw                         Raw data output from minion
   basefolder/2_SampleList_and_Rampart      List of all samples
   basefolder/3_Artic_Output                Processed files from artic pipeline
   basefolder/4_Consensus                   Consensus sequences multi-fasta per run
   basefolder/5_GISAID                      Sequences for uploading to GISAID
   
   Syntax: $(basename "$0") [-123|b|n|p|r|]
      
    Mandatory args:
    -b      Basefolder for all inputs and outputs - please use an absolute path
    -n      Runname to prefix output files, folder names etc

    Optional args:
    -h      Print this help
    -m      Process using the medaka pipeline rather than with nanopolish
    -p      Primers used during PCR [Sanger | Artic]
    -r      Location directory of the primer scheme used

    Optional steps to omit:
    -1      Omit Step 1 (guppy_barcoder)
    -2      Omit Step 2 (guppyplex)
    -3      Omit Step 3 (artic minion)
    -4      Omit Step 4 (Consensus sequences)
    -5      Omit Step 5 (Sequencing statistics)
   "
}

function present {
    if [ $2 = 'd' ] ; then
        if [ ! -d $1 ] ; then
            printf "${RED}ERROR:${NC} $1 is not present \nExiting script."
            exit
        fi
    elif [ $2 = 'f' ] ; then
        if [ ! -f $1 ] ; then
            printf "${RED}ERROR:${NC} $1 is not present \nExiting script."
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
    if [ $1 = 1 ] ; then
        Help
        printf "${RED}ERROR:${NC} $2 variable not supplied. See Help message above.\n Exiting script\n\n"
        exit   
    fi
}

#==============Recognise any CLI options==================
while getopts 'hp:b:r:n:m12345' opt; do
  case "$opt" in
    b)
      BASEFOLDER=$OPTARG
      ;;
    p)
      PRIMERS=$OPTARG
      ;;
    r)
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
    n)
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
#ANSI escape codes: 
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
RED='\033[0;31m'
GREEN='\033[0;32m' 
BLUE='\033[0;34m' 
ORANGE='\033[0;33m' 
NC='\033[0m' # No Color

#Create a space between cmd
printf "\n${ORANGE}##### Running ArticProcess.sh ##### ${NC}\n"

#################################################
#Enter default values if parameter not defined or build from exisiting args
# These are MANDATORY
BASEFOLDER=${BASEFOLDER:- 1}
check_var $BASEFOLDER "-b (basefolder)"
BASEFOLDER=`realpath $BASEFOLDER`

present $BASEFOLDER "d"

RUNNAME=${RUNNAME:- 1}
check_var $RUNNAME "-n (runname)"

S1=${S1:- 1}
S2=${S2:- 1}
S3=${S3:- 1}
S4=${S4:- 1}
S5=${S5:- 1}

#################################################
#Check artic environment present
if [ $(which artic | wc -l) = 0 ] ; then
    printf "${RED}ERROR:${NC} artic environment not present. Please activate the appropriate environment e.g.:\n\n    
    conda activate artic\n    
    \nExiting script\n"
    exit 
fi

#################################################
#Determine if the necessary files and directories and arguments etc exist or create them if not


#Make directories for the demultiplexed, combined and size selected fast_q files
ARTIC_OUT="$BASEFOLDER/3_Artic_Output/$RUNNAME"
check_mkdir $ARTIC_OUT
check_mkdir "$ARTIC_OUT/fastq"
check_mkdir "$ARTIC_OUT/processed"
RAWDATADIR="$BASEFOLDER/1_Raw/$RUNNAME"
present $RAWDATADIR "d"
SEQUENCINGSUMMARY=$(find $RAWDATADIR -name 'sequencing_summary*' -not -name '*.tmp')
if [ -z "$SEQUENCINGSUMMARY" ] ; then
    printf "${RED}ERROR:${NC} No sequencing summary found for $RUNNAME run.\n"
fi

###############
#Guppy Barcoder (step 1)
if [ $S1 = 1 ] ; then
    #Check guppy_barcoder present
    if [ $(which guppy_barcoder | wc -l) = 0 ] ; then
        printf "${RED}ERROR:${NC} guppy_barcoder not present. Exiting script\n\n"
        exit
    fi
    FASTQRAW="$RAWDATADIR/fastq_pass"
    present $FASTQRAW "d"
fi

###############
#Guppy Plex (step 2) and step 3 identical
if [ $S2 = 1 ] || [ $S3 = 1 ]; then
    #Check primers and determine max and min sizes and primer scheme for analysis
    PRIMERS=${PRIMERS:- 1}
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
    MEDAKA=${MEDAKA:- 0}
    
    #Check primer scheme supplied
    PRIMERDIR=${PRIMERDIR:- 1}
    check_var $PRIMERDIR "-r (primer scheme location)"
    PRIMERDIR=`realpath $PRIMERDIR`
    
    #Identify csv file
    SAMPLEFILE="$BASEFOLDER/2_SampleList_and_Rampart/Samples_Sequenced.csv"
    present $SAMPLEFILE "f"
    #Read in Samples
    readarray -t SAMPLES < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$2}' $SAMPLEFILE)
    #Check for sample duplicates
    if [ `printf '%s\n' "${SAMPLES[@]}" | awk '!($0 in seen){seen[$0];next} 1' | wc -l` -gt 0 ] ; then
        printf "${RED}ERROR:${NC} The following duplicate samples are present in $SAMPLEFILE file\n"
        printf '%s\n' "${SAMPLES[@]}" | awk '!($0 in seen){seen[$0];next} 1'
        exit
    fi
    #Read in barcodes
    readarray -t BARCODES < <(awk -F"," '$7 ~ /'$RUNNAME'/ {print$1}' $SAMPLEFILE)
    #Check for barcode duplicates
    if [ `printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1' | wc -l` -gt 0 ] ; then
        printf "${RED}ERROR:${NC} The following duplicate barcodes are present in $SAMPLEFILE file\n"
        printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1'
        exit
    fi
fi

printf "${GREEN}CHECKED:${NC}All required programs, files and locations are present.
${GREEN}CHECKED:${NC}No errors in sample list.\n"

#==============THE SCRIPT==================Se   
#STEP 1: Run the guppy barcoder to demultiplex into separate barcodes
if [ $S1 = 1 ] ; then
    printf "\n###### ${BLUE}Step 1: Running the guppy_barcoder to demultiplex the FASTQ files.${NC} ######\n\n"
    guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
    printf "\n###### ${GREEN}Step 1: guppy_barcoder completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 1: Skipping guppy_barcoder step${NC} ######\n\n"
fi

#STEP 2: Combine all identical barcode reads into same fasta and exclude sizes
if [ $S2 = 1 ] ; then
    printf "\n###### ${BLUE}Step 2: Combining demultiplexed files into a single fastq and excluding based on size.${NC} ######\n\n"
    readarray -t S2DIRS < <(find $ARTIC_OUT/fastq -type d -name 'barcode[0-9]*')
    
    #Change directory as unable to redirect output from guppyplex
    cd $ARTIC_OUT/fastq
    
    #Run through the demuxed barcodes to combine into single fastq files
    for S2DIR in "${S2DIRS[@]}" ; do
       artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME
    done
    echo -e "\n###### ${GREEN}Step 2: guppyplex completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 2: Skipping guppyplex step${NC} ######\n\n"
fi

#STEP 3: Processing with artic minion pipeline
if [ $S3 = 1 ] ; then
    printf "\n###### ${BLUE}Step 3: Importing samplenames and processing with artic minion command.${NC} ######\n\n"
    
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
            printf "${RED}WARNING:${NC} $FILE is not present - moving to next sample.\n"
        else
            #Ensure fastq file is not empty
            FASTQLENGTH=`wc -l $FILE | sed s/$FILE//`  
            
            #Remove files that have not got enough data
            if [ $FASTQLENGTH -gt 1000 ] ; then
                echo -e "\n\n \n${ORANGE} Processing barcode number $BARCODE, Sample $SAMPLENAME from file $FILE \n${NC}"
                #Run processing scheme
                if [ $MEDAKA = 1 ] ; then
                    printf "\n${GREEN}Using Medaka pipeline${NC}\n"
                    artic minion --medaka --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE $PRIMERSCHEME $SAMPLENAME
                else
                    printf "\n${GREEN}Using Nanopolish pipeline${NC}\n"
                    artic minion --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME
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
                printf "${RED}ERROR:${NC} Too few reads (n = $FASTQLENGTH) in File $FILE (Barcode $BARCODE, Sample $SAMPLENAME).
                Aborting processing this file\n"
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
    find /$ARTIC_OUT/processed -name '*.consensus.fasta' | xargs cat > $RUNNAME.consensus.fasta
    
    echo -e "\n###### ${GREEN}Step 4: consensus sequences compiled. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 4: Skipping consensus sequences ${NC} ######\n\n"
fi

#STEP 5: Generate GISAID stats
if [ $S5 = 1 ] ; then
    printf "\n###### ${BLUE}Step 5: Generating statistics on sequencing. ${NC} ######\n\n"
    #Run 
    python run_gisaid-statistics.py -b /home/dan/TESTWGS
    echo -e "\n###### ${GREEN}Step 5: Sequencing statistics compiled. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 5: Skipping generating sequencing statistics${NC} ######\n\n"
fi


exit

######################## NOTES ##########

# Concatenate and tidy up the names of the fasta headers
cat *.fasta | sed 's/\/ARTIC\/nanopolish MN908947.3/\/2020/'| sed 's/>/>SARS-CoV-2\/human\/Zambia\//' > out.fasta

#Bug catching - list all defined variables
( set -o posix ; set ) | less
