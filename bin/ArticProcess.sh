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

VERSION="1.1"
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
   basefolder/6_QCAnalysis                  QC output from ncov-tools
   
   Syntax: $(basename "$0") [-123456hv|b|r|m|p|s]
   Version: $VERSION
   
    Mandatory args:
    -b      Basefolder for all inputs and outputs - please use an absolute path
    -r      Runname to prefix output files, folder names etc

    Optional args:
    -h      Print this help
    -m      Process using the medaka pipeline rather than with nanopolish
    -p      Primers used during PCR [Sanger | Artic | Midnight]
    -s      Location directory for github clones (assumes Path-nCoV and ncov-tools cloned into this location)
    -v      Version

    Optional steps to omit:
    -1      Omit Step 1 (guppy_barcoder)
    -2      Omit Step 2 (guppyplex)
    -3      Omit Step 3 (artic minion)
    -4      Omit Step 4 (Consensus sequences)
    -5      Omit Step 5 (Sequencing statistics)
    -6      Omit Step 6 (Generate JSON file for RAMPART)
    -7      Omit Step 7 (ncov-tools QC pipeline)
    "
    
    HelpInstall
}

function HelpInstall {
    printf "\n${RED}   INSTALLATION:${NC}
    For all parts of the script to run, the following environments need to be installed:
        ${GREEN}artic${NC} - https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
        ${GREEN}pangolin${NC} - https://cov-lineages.org/resources/pangolin/installation.html
        ${GREEN}nextclade${NC} - https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html
        ${GREEN}ncov-tools${NC} - https://github.com/jts/ncov-tools
    
    Other tools that need to be in the users path are:
        ${GREEN}sequencing_statistics.py${NC} - from the Path_nCoV git repository
        ${GREEN}guppy${NC} - https://community.nanoporetech.com/downloads (requires ONT community login)
        
    Finally for nextclade, the datasets are expected to be in:
        ~/.nextclade/dataset/sars-cov-2/
    OR the git repo (-s flag) location (may not be most up-to-date)
    (see https://docs.nextstrain.org/projects/nextclade/en/latest/user/datasets.html for more info)
    jq should be installed to allow the script to check dataset and nextcladeCli compatibility\n\n"
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
    printf "   ${BLUE}$1${NC} found\n"
}

function check_mkdir {
    if [ ! -d $1 ] ; then
        echo -e "Making directory $1"
        mkdir -p $1
    fi
}

function check_var {
#$1 variable required
#$2 Help message identifying what variable relates to
    if [ -z $1 ] || [ $1 = 1 ]; then
        Help
        printf "\n${RED}ERROR:${NC} $2 variable not supplied. See Help message above.\n Exiting script\n\n"
        exit   
    fi
}

function check_package {
    if [ $(which $1 | wc -l) = 0 ] ; then
        printf "   ${RED}ERROR:${GREEN} $1${NC} executable not found.\n\n"
        HelpInstall
        printf "${RED}Exiting script${NC}\n\n"
        exit
    else
        printf "      ${GREEN}$1${NC} command found\n"
    fi
}

function GETCOLMID {
#Function to identify the column number for a header
awk -F"$1" 'NR==1{ for (i=1;i<=NF;i++) if ($i == "'$2'") print i }' $3
}

function check_condaenv {
#Test if environment is present
readarray -t LIST < <(conda info --envs | grep envs | awk '{print $1}' | grep -e ^$1)
if [ ${#LIST[@]} = 0 ] ; then 
    printf "   ${RED}ERROR:${GREEN} $1${NC} conda environment not found\n"
    exit
elif [ ${#LIST[@]} = 1 ] ; then #Use if single match
    printf "   ${GREEN}$1${NC} environment found\n"
#    ARTIC=$1
elif [ ${#LIST[@]} > 1 ] ; then
    printf "   Found multiple matches for $1:\n"
    for L in "${LIST[@]}" ; do
        printf "     ${GREEN}$L${NC}\n"
    done
fi
}

function change_conda {
    CONDA=`which conda | sed s'/\/condabin//' | sed s'/bin\///' | sed s'/\/conda//'`
    source $CONDA/etc/profile.d/conda.sh 
    set +eu
    conda activate $1
}

#==============Recognise any CLI options==================
while getopts 'vhp:b:r:s:m1234567' opt; do
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
      GITDIR=$OPTARG
      ;;
    h)
      Help
      exit
      ;;
    m)
      MEDAKA=1
      MEDAKAMODEL="r941_min_high_g303"
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
    7)
      S7=0
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
printf "\n${GREEN}##### Running `basename $0` Version $VERSION ##### ${NC}\n"
printf "${GREEN}Checking all files, folders, environments and data are present${NC}\n"
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
S7=${S7:- 1}
MEDAKA=${MEDAKA:- 0}            # 3
REFDIR=${REFDIR:- 1}      # 3
PRIMERS=${PRIMERS:- 1}          # 2, 3
ARTIC_OUT="$BASEFOLDER/3_Artic_Output/$RUNNAME"
RAWDATADIR="$BASEFOLDER/1_Raw/$RUNNAME"
GITDIR=${GITDIR:- 1}

#################################################
#Make directories for the demultiplexed, combined and size selected fast_q files

#These are folders that are absolutely required
DIRS_REQ=("$RAWDATADIR" "$BASEFOLDER/2_SampleList_and_Rampart")
for DIR_REQ in "${DIRS_REQ[@]}" ; do
    present $DIR_REQ "d"
done

#These folders may only be created during the run
DIRS_ADD=("$ARTIC_OUT/fastq" "$ARTIC_OUT/processed" "$BASEFOLDER/4_Consensus" "$BASEFOLDER/5_GISAID" "$BASEFOLDER/5_GISAID/intermediates" "$BASEFOLDER/5_GISAID/proteins" "$BASEFOLDER/6_QCAnalysis") 
for DIR_ADD in "${DIRS_ADD[@]}" ; do
    check_mkdir $DIR_ADD
done

ALLSEQ="allsequences.fasta" #Name of file for all sequences to output to

#Determine if the necessary files and directories and arguments etc exist for each step
printf ""
if [ $S1 = 1 ] ; then
    check_package "guppy_barcoder"
    FASTQRAW="$RAWDATADIR/fastq_pass"
    present $FASTQRAW "d"
fi

if [ $S2 = 1 ] || [ $S3 = 1 ]; then
    check_condaenv artic
    change_conda artic
    check_package artic
fi

if [ $S2 = 1 ] || [ $S3 = 1 ] || [ $S7 = 1 ]; then
    #Check primers and determine max and min sizes and primer scheme for analysis
    check_var $PRIMERS "-p (primers)"
    if [ $PRIMERS = "Sanger" ] ; then
        MIN=750
        MAX=1250
        PRIMERSCHEME=SARS-CoV-2/V1
    elif [ $PRIMERS = "Artic" ] ; then
        MIN=400
        MAX=700
        PRIMERSCHEME=SARS-CoV-2/V3
    elif [ $PRIMERS = "Midnight" ] ; then
        MIN=950
        MAX=1450
        PRIMERSCHEME=SARS-CoV-2/V2
    else
        Help
        printf "${RED}ERROR:${NC} Unrecognised primer scheme\n\n"
        exit
    fi
    printf "   Primerscheme - $PRIMERSCHEME (Min: $MIN, Max: $MAX)\n"
fi

if [ $S3 = 1 ] || [ $S5 = 1 ] || [ $S7 = 1 ]; then
    #Check primer scheme supplied
    check_var $GITDIR "-s (Git projects directory)"
    REFDIR=`realpath "$GITDIR/Path_nCoV/reference"`
    NCTDIR=`realpath "$GITDIR/ncov-tools"`
    present "$REFDIR" "d"
    present "$NCTDIR" "d"
fi

if [ $S5 = 1 ] ; then
    #Pangolin requirements
    check_condaenv pangolin
    change_conda pangolin
    check_package pangolin
    
    #nextclade requirements
    check_condaenv nextclade
    change_conda nextclade
    check_package nextclade
    
    #Identify dataset location to use
    NCLADE_DATA=`echo $HOME/.nextclade/dataset/sars-cov-2/`
    if [ ! -d $NCLADE_DATA ] ; then
      NCLADE_DATA="$REFDIR/nextclade/dataset"
    fi
    present $NCLADE_DATA "d"
    #Check whether dataset matches installed nextclade version
    if [ $(which jq | wc -l) = 1 ] ; then
        printf "   ${BLUE}Identified jq - checking nextclade compatibility${NC}\n"
        NC_DATASET=`jq '.compatibility.nextcladeCli.min' $NCLADE_DATA/tag.json | sed s/\"//g`
        NC_VERSION=`nextclade -v`
        if [ $(printf "${NC_VERSION}\n${NC_DATASET}" | sort -V | head -1) = "${NC_VERSION}" ]; then
            printf "${ORANGE}WARNING:${NC} Nextclade dataset requires minimum v$NC_DATASET nextcladeCli ($NC_VERSION installed). Upgrade strongly suggested\n"
        else
            printf "   ${BLUE}Nextclade version $NC_VERSION used with dataset requiring $NC_DATASET\n${NC}"
        fi
    else
        printf "   ${ORANGE}WARNING:${NC} jq not installed unable to check compatibility of nextclade dataset and installed version${NC}\n"
    fi    
     
    #Sequencing stats requirements
    check_package sequencing_statistics.py
fi

if [ $S7 = 1 ] ; then
    check_condaenv ncov-qc
    change_conda ncov-qc
    check_package snakemake
fi

if [ $S3 = 1 ] ; then
    SEQUENCINGSUMMARY=$(find $RAWDATADIR -name 'sequencing_summary*' -not -name '*.tmp')
    check_var ${SEQUENCINGSUMMARY:- 1} "Sequencing Summary"
    printf "   Sequencing summary file present\n"
fi

if [ $S3 = 1 ] || [ $S6 = 1 ] || [ $S7 = 1 ] ; then
    #Identify csv file
    SAMPLEFILE="$BASEFOLDER/2_SampleList_and_Rampart/Samples_Sequenced.csv"
    present $SAMPLEFILE "f"
    printf "   Samplefile present\n"
fi


if [ $S3 = 1 ] || [ $S6 = 1 ] ; then
    #Check if all of the seqIDs are unique in the samples file
    #Identify correct column IDs
    SEQIDCOL="$(GETCOLMID "," "SeqID" "$SAMPLEFILE")"
    SEQRUNCOL="$(GETCOLMID "," "SeqRun" "$SAMPLEFILE")"
    SEQBARCODE="$(GETCOLMID "," "SeqBarcode" "$SAMPLEFILE")"
    
    readarray -t IDS < <(awk -F"," -v i=$SEQIDCOL 'NR>1 {if ($i) print$i}' $SAMPLEFILE | sort | uniq -d)
    if [ ${#IDS[@]} -gt 0 ]; then
        printf "${RED}ERROR:${NC} The following duplicate seqIDs exist in $SAMPLEFILE. 
        Please correct so that every seqID (column $SEQIDCOL) is unique to the entire list (empty records are not included).\n"
        printf '%s\n' "${IDS[@]}"
        exit
    else
        printf "   Sequence ID array accepted\n"
    fi
  
    #Read in Samples
    readarray -t SAMPLES < <(awk -F"," -v i=$SEQRUNCOL -v j=$SEQIDCOL '$i ~ /'$RUNNAME'/ {print$j}' $SAMPLEFILE)
    if (( ${#SAMPLES[@]} == 0 )); then
        printf "\n${RED}ERROR:${NC} SAMPLES array is empty for run $RUNNAME\nExiting script\n\n"
        exit
    else
        printf "   Sample ID array accepted\n"
    fi
    
    #Read in barcodes
    #readarray -t BARCODES < <(awk -F"," '$10 ~ /'$RUNNAME'/ {print$8}' $SAMPLEFILE)
    readarray -t BARCODES < <(awk -F"," -v i=$SEQRUNCOL -v j=$SEQBARCODE '$i ~ /'$RUNNAME'/ {print$j}' $SAMPLEFILE)
    if (( ${#BARCODES[@]} == 0 )); then
        printf "\n${RED}ERROR:${NC} BARCODES array is empty \nExiting script\n\n"
        exit
    fi
    
    #Check for barcode duplicates within the run
    if [ `printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1' | wc -l` -gt 0 ] ; then
        printf "${RED}ERROR:${NC} The following duplicate barcodes are present in $SAMPLEFILE file for $RUNNAME run:\n"
        printf '%s\n' "${BARCODES[@]}" | awk '!($0 in seen){seen[$0];next} 1'
        exit
    else
        printf "   Sample barcode array accepted\n"
    fi
    printf "   ${GREEN}Sample List:${NC} No errors identified.\n"
fi

#Define logfile output
LOGFOLDER="$BASEFOLDER/Logs/"
RUNLOG="$LOGFOLDER/$RUNNAME/"
check_mkdir $RUNLOG
RUNLOG=$(echo "$RUNLOG${RUNNAME}_Step_")

#All checks complete
printf "${GREEN}CHECKED:${NC}All required programs, files and locations are present.\n\n"

#==============THE SCRIPT==================Se   
#STEP 1: Run the guppy barcoder to demultiplex into separate barcodes
if [ $S1 = 1 ] ; then
    printf "\n###### ${BLUE}Step 1: Running the guppy_barcoder to demultiplex the FASTQ files.${NC} ######\n\n"
    
    if [ `find $FASTQRAW -type d -name 'barcode*' | wc -l ` -gt 0 ] ; then
        printf "Found barcoded data structure\n"
        if [ `find $FASTQRAW -type f -name '*.fastq.gz' | wc -l ` -gt 0 ] ; then
            printf "Unzipping files saved as a gzip\n"
            find $FASTQRAW -type f -name '*.fastq.gz' | xargs gunzip
        fi
        printf "Moving all fastq files to $FASTQRAW directory and deleting subdirectories\n"
        find $FASTQRAW -type f -name '*.fastq' | xargs -I '{}'  mv {} ${FASTQRAW}/
        rm -rvf $FASTQRAW/barcode*
        rm -rvf $FASTQRAW/unclassified
    else
        printf "Found standard data structure\n"
    fi
    
    #Test for version of guppy to determine command to run
    GUPVER=`guppy_barcoder -v | grep Version | cut -c 101-101`
    printf "Identified guppy_barcoder version $GUPVER\n"
    
    if [ $GUPVER -lt 6 ] ; then
        printf "guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files barcode_arrs_nb96.cfg${NC}\n" | tee "${RUNLOG}1.log"
        guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files "barcode_arrs_nb96.cfg" | tee -a "${RUNLOG}1.log"
    else
        printf "guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --barcode_kits SQK-NBD112-96 ${NC}\n" | tee "${RUNLOG}1.log"
        guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --barcode_kits "SQK-NBD112-96"  | tee -a "${RUNLOG}1.log"
    fi
    
    printf "\n###### ${GREEN}Step 1: guppy_barcoder completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 1: Skipping guppy_barcoder step${NC} ######\n\n"
fi

#STEP 2: Combine all identical barcode reads into same fasta and exclude sizes
if [ $S2 = 1 ] ; then
    printf "\n###### ${BLUE}Step 2: Combining demultiplexed files into a single fastq and excluding based on size.${NC} ######\n\n" | tee "${RUNLOG}2.log"
    readarray -t S2DIRS < <(find $ARTIC_OUT/fastq -type d -name 'barcode[0-9]*')
    
    change_conda artic
    #Change directory as unable to redirect output from guppyplex
    cd $ARTIC_OUT/fastq
    
    #Run through the demuxed barcodes to combine into single fastq files
    for S2DIR in "${S2DIRS[@]}" ; do
        printf "${LG}artic guppyplex --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME ${NC}\n" | tee -a "${RUNLOG}2.log" #--skip-quality-check
       artic guppyplex --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME | tee -a "${RUNLOG}2.log"
    done
    echo -e "\n###### ${GREEN}Step 2: guppyplex completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 2: Skipping guppyplex step${NC} ######\n\n"
fi

#STEP 3: Processing with artic minion pipeline
if [ $S3 = 1 ] ; then
    printf "\n###### ${BLUE}Step 3: Importing samplenames and processing with artic minion command.${NC} ######\n\n" | tee "${RUNLOG}3.log"
    
    change_conda artic
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
            printf "${RED}WARNING:${NC} $FILE is not present - moving to next sample.\n" | tee -a "${RUNLOG}3.log"
        else
            #Ensure fastq file is not empty
            FASTQLENGTH=`wc -l $FILE | sed s/$FILE//`  
            
            echo -e "\n\n \n${ORANGE} Processing barcode number $BARCODE, Sample $SAMPLENAME from file $FILE \n${NC}" | tee -a "${RUNLOG}3.log"
            
            #Remove files that have not got enough data
            if [ $FASTQLENGTH -gt 1 ] ; then
                #Run processing scheme
                if [ $MEDAKA = 1 ] ; then
                    printf "\n${GREEN}Using Medaka pipeline, with the Medaka Model $MEDAKAMODEL ${NC}\n"
                    printf "${LG}artic minion --medaka --medaka-model $MEDAKAMODEL --normalise 200 --threads 24 --scheme-directory $REFDIR/primer-schemes/ --read-file $FILE $PRIMERSCHEME $SAMPLENAME ${NC}\n" | tee -a "${RUNLOG}3.log"
                    artic minion --medaka --medaka-model $MEDAKAMODEL --normalise 200 --threads 24 --scheme-directory $REFDIR/primer-schemes/ --read-file $FILE $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${RUNLOG}3.log"
                else
                    printf "\n${GREEN}Using Nanopolish pipeline${NC}\n"
                    printf "${LG}artic minion --normalise 200 --threads 24 --scheme-directory $REFDIR/primer-schemes/ --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME${NC}\n" | tee -a "${RUNLOG}3.log"
                    artic minion --normalise 200 --threads 24 --scheme-directory $REFDIR/primer-schemes/ --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${RUNLOG}3.log"
                fi
                
                #Make samplename subdirectory if required
                check_mkdir $ARTIC_OUT/processed/$SAMPLENAME
                #Move output from what was produced into its own directory in processed
                find ./ -type f -name "$SAMPLENAME*" | xargs -I '{}'  mv {} "$ARTIC_OUT/processed/$SAMPLENAME"/
            else
                printf "\n\n${RED}ERROR:${NC} Too few reads (n = $FASTQLENGTH) in File $FILE (Barcode $BARCODE, Sample $SAMPLENAME).\nAborting processing this file\n" | tee -a "${RUNLOG}3.log"
            fi
        fi
    done
    echo -e "\n###### ${GREEN}Step 3: artic minion completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 3: Skipping artic minion step${NC} ######\n\n"
fi

#STEP 4: Concatenate consensus sequences
if [ $S4 = 1 ] ; then
    printf "\n###### ${BLUE}Step 4: Concatenating consensus sequences from artic pipeline${NC} ######\n\n"
    cd $BASEFOLDER
    #Ensure folder exists
    CONSENSUS="$BASEFOLDER/4_Consensus"
    check_mkdir $CONSENSUS
    cd $CONSENSUS
    
    #Concatenate all of the fasta consensus data from this run
    find /$ARTIC_OUT/processed -name '*.consensus.fasta' | xargs cat > $RUNNAME.consensus.fasta
    #Cat all runs sequences together
    find ./ -type f -name "*consensus.fasta" | xargs -I '{}'  cat {} > "sequences.fasta"
    #Remove the additional info in the header to just leave the Sample ID
    cat sequences.fasta | sed  s'/\/ARTIC\/nanopolish MN908947.3//' > $ALLSEQ
    rm sequences.fasta
    
    echo -e "\n###### ${GREEN}Step 4: Consensus sequences compiled ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 4: Skipping consensus sequences ${NC} ######\n\n"
fi

#STEP 5: Generate stats
if [ $S5 = 1 ] ; then
    printf "\n###### ${BLUE}Step 5: Generating statistics and lineages on sequences. ${NC} ######\n"
    printf "\n###### ${GREEN} Copying all consensus sequences (${ALLSEQ}) ${NC} ######\n\n"
    cd "$BASEFOLDER/5_GISAID"
    cp ../4_Consensus/$ALLSEQ ./
    printf "Done"
    
    printf "\n###### ${GREEN} Determining nextclade lineages ${NC} ######\n\n"
    change_conda nextclade
    
    printf "nextclade --verbose --in-order --input-fasta $ALLSEQ --input-dataset $NCLADE_DATA --input-pcr-primers $REFDIR/nextclade/primers_UNZA.csv --output-csv nextclade.csv --output-tree nextclade.auspice.json --output-basename allsequences \n\n" | tee "${LOGFOLDER}nextclade.log"
    nextclade --verbose --in-order --input-fasta $ALLSEQ --input-dataset $NCLADE_DATA --input-pcr-primers $REFDIR/nextclade/primers_UNZA.csv --output-csv nextclade.csv --output-tree nextclade.auspice.json --output-basename allsequences 2>&1 | tee -a "${LOGFOLDER}nextclade.log"
    #Check everything ran properly
    if [ `grep -c ERROR ${LOGFOLDER}nextclade.log` != 0 ] ; then
        printf "${ORANGE}WARNING: Error(s) identified in Nextclade log\n"
    fi
    #Clean-up nextclade outputs
    mv allsequences.gene* proteins/ 
    mv nextclade* intermediates/
    mv allsequences*.csv intermediates/

    printf "\n###### ${GREEN} Determining PANGO lineages ${NC} ######\n\n"
    change_conda pangolin
    pangolin $ALLSEQ 2>&1 | tee "${LOGFOLDER}pango.log"
    mv lineage_report.csv intermediates/
        
    printf "\n###### ${GREEN} Determining sequencing statistics, merging with PANGO, Nextclade and metadata ${NC} ######\n\n"
    sequencing_statistics.py -d $BASEFOLDER 2>&1 | tee "${LOGFOLDER}gisaid.log"
        
    printf "\n###### ${GREEN}Step 5: Sequencing statistics compiled. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 5: Skipping generating sequencing statistics${NC} ######\n\n"
fi

if [ $S6 = 1 ] ; then
    printf "\n###### ${BLUE}Step 6: Generating JSON files for Rampart. ${NC} ######\n\n"

    COUNT=0
    
    #Store output in run specific folder
    RAMPARTOUTPUT="$BASEFOLDER/2_SampleList_and_Rampart/${RUNNAME}/"
    check_mkdir $RAMPARTOUTPUT
    JSONFILE="${RAMPARTOUTPUT}run_configuration.json"
    #Start the JSON File
    printf "{
    \"title\": \"Run $RUNNAME\",
    \"basecalledPath\": \"fastq_pass\",
    \"samples\": [" > $JSONFILE
    
    #Rampart can only deal with the first 24 barcodes natively. If >24 then rewrite output to recognise guppy output
    if [ ${#BARCODES[@]} -gt 24 ] ; then
        PREFIX="barcode"
    else
        PREFIX="NB"
    fi

    for BARCODE in "${BARCODES[@]}"; do 
        SAMPLE=${SAMPLES[$COUNT]}
        #Ensure that all barcodes are 2 digit
        if [ ${#BARCODE} -lt 2 ]; then
            BARCODE="0$BARCODE"
        fi
        printf "    {
        \"name\": \"$SAMPLE ($BARCODE)\",
        \"description\": \"\",
        \"barcodes\": [ \"${PREFIX}${BARCODE}\" ]" >> $JSONFILE
        if [ ${#BARCODES[@]} == $((COUNT+1)) ]; then
            echo -e "    }" >> $JSONFILE #If last entry then there shouldn't be a comma
        else
            echo -e "    }," >> $JSONFILE 
        fi
        #Increment count
        (( COUNT += 1 )) 
        printf "${BLUE}ADDED:${NC} Sample = $SAMPLE, Barcode = $BARCODE\n"
    done
    #Finish off the JSON File
    printf "  ]\n}" >> $JSONFILE
    printf "\n###### ${GREEN}Step 6: Rampart config saved to $JSONFILE. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 6: Skipping Generating JSON files for Rampart${NC} ######\n\n"
fi

if [ $S7 = 1 ] ; then
    printf "###### ${BLUE}Step 7: Running ncov-tools QC pipeline${NC} ######\n\n"
    change_conda "ncov-qc"
    QCDIR="$BASEFOLDER/6_QCAnalysis/$RUNNAME"
    check_mkdir $QCDIR
       
    cd "$BASEFOLDER/6_QCAnalysis/$RUNNAME"
    
    printf "Removing prior QC run details:\n"
    rm -rvf $QCDIR/*
    
    printf "Generating Metadata.tsv file\n"
    METADATAFN="Metadata.tsv"
    SEQIDCOL="$(GETCOLMID "," "SeqID" "$SAMPLEFILE")"
    CTCOL="$(GETCOLMID "," "Ctvalue" "$SAMPLEFILE")"
    SEQRUNCOL="$(GETCOLMID "," "SeqRun" "$SAMPLEFILE")"
    awk -F"," -v i=$SEQRUNCOL -v k=$SEQIDCOL -v l=$CTCOL 'NR==1 {print "sample","ct"} $i ~ /'$RUNNAME'/ {print$k,$l}' OFS='\t' $SAMPLEFILE > $METADATAFN
    NTCLIST=`awk -F"\t" '{print$1}' $METADATAFN | grep NTC | sed -e 's/^/\\"/' | sed -e 's/$/\\"/' | awk -vORS=, '{print$1}'`

    printf "\nGenerating config.yaml file\n"
    printf "data_root: $BASEFOLDER/3_Artic_Output/$RUNNAME/processed
run_name: $RUNNAME
reference_genome: ${REFDIR}/primer-schemes/${PRIMERSCHEME}/SARS-CoV-2.reference.fasta
platform: oxford-nanopore
primer_bed: ${REFDIR}/primer-schemes/${PRIMERSCHEME}/SARS-CoV-2.bed
bam_pattern: \"{data_root}/{sample}/{sample}.sorted.bam\"
consensus_pattern: \"{data_root}/{sample}/{sample}.consensus.fasta\"
variants_pattern: \"{data_root}/{sample}/{sample}.pass.vcf.gz\"
skip_empty_negatives: true
metadata: \"$METADATAFN\"
negative_control_samples: [ $NTCLIST ]\n" > "config.yaml"

    printf "Running snakemake report\n"
    snakemake --cores all -s $NCTDIR/workflow/Snakefile all_final_report
    
else
    printf "###### ${GREEN}Step 7: Skipping generation of QC stats${NC} ######\n\n"
fi
exit
