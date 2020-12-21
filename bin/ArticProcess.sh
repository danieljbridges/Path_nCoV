#!/bin/bash
# ----------------------------------------------------------------------
# Script to automate the artic processing pipeline for SARS-CoV-2 WGS
# Version 0.2 29/10/2020
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
   echo -e "This script chains together the various processes of the artic network bioinformatic pipeline for each sample.
   Syntax: $(basename "$0") [-f123|n|o|p|r|s|w]
      
    Configuration:
   -n      Runname to prefix output files with (will be concatenated with the sample name)
   -o      Directory location where the output files should be stored
   -p      Primers used during PCR [Sanger | Artic]
   -r      Location directory of the primer scheme used
   -w      Directory location of the raw output from the ONT platform (assumed to contain fast5_pass and fastq_pass directories.
   
    Optional:
   -h      Print this help
   -s      Name of tab delimited file containing the barcode (2 digit) and sample name relationship, if not Samples.txt
   -m      Process using the medaka pipeline rather than with nanopolish
   -1      Omit the guppy_barcoder step
   -2      Omit the guppyplex step
   -3      Omit the artic minion step
   "
}

function present {
    if [ $2 = 'd' ] ; then
        if [ ! -d $1 ] ; then
            echo -e "$1 is not present \nExiting script."
            exit
        fi
    elif [ $2 = 'f' ] ; then
        if [ ! -f $1 ] ; then
            echo -e "$1 is not present \nExiting script."
            exit
        fi
    fi
}

#==============Recognise any CLI options==================
while getopts 'hp:o:r:n:s:w:123' opt; do
  case "$opt" in
    p)
      PRIMERS=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    r)
      PRIMERDIR=$OPTARG
      present $PRIMERDIR "d"
      ;;
    w)
      RAWDATADIR=$OPTARG
      FAST5RAW="$RAWDATADIR" #/fast5_pass"
      FASTQRAW="$RAWDATADIR/fastq_pass"
      #Ensure all present
      present $FAST5RAW "d"
      present $FASTQRAW "d"
      present $RAWDATADIR "d"
      ;;
    h)
      Help
      exit
      ;;
    s)
      SAMPLES=$OPTARG
      present $SAMPLES "f"
      ;;
    m)
      MEDAKA=1
      ;;
    1)
      GUPPYBARCODER=0
      ;;
    2)
      GUPPYPLEX=0
      ;;
    3)
      ARTICMINION=0
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

#Check artic environment present
if [ $(which artic | wc -l) = 0 ] ; then
    printf "${RED}ERROR:${NC} artic environment not present. Please activate the appropriate environment:\n\n    
    conda activate artic-ncov2019\n    
    \nExiting script\n"
    exit 
fi
#Check guppy_barcoder present
if [ $(which guppy_barcoder | wc -l) = 0 ] ; then
    printf "${RED}ERROR:${NC} guppy_barcoder not present. Exiting script\n\n"
    exit
fi

#Determine if the files, directories etc exist 
#Find the name of the sequencing summary
SEQUENCINGSUMMARY="$RAWDATADIR/`ls $RAWDATADIR | grep sequencing_summary`"
present $SEQUENCINGSUMMARY "f"

#Enter default values if parameter not defined
SAMPLES=${SAMPLES:-$RAWDATADIR/Samples.txt}
GUPPYBARCODER=${GUPPYBARCODER:- 1}
GUPPYPLEX=${GUPPYPLEX:- 1}
ARTICMINION=${ARTICMINION:- 1}
MEDAKA=${MEDAKA:- 0}

#Make output dir if needed
if [ ! -d $OUTPUT ] ; then
    mkdir $OUTPUT
    echo "Making $OUTPUT directory"
fi

#Make a subdirectory for the demultiplexed, combined and size selected fast_q files
FASTQ_OUT="${OUTPUT}/fastq"
if [ ! -d $FASTQ_OUT ] ; then
    echo -e "Making directory $FASTQ_OUT"
    mkdir $FASTQ_OUT
fi

#Make a subdirectory for the processed output from the pipeline
FINAL_OUT="${OUTPUT}/processed"
if [ ! -d $FINAL_OUT ] ; then
    echo -e "Making directory $FINAL_OUT"
    mkdir $FINAL_OUT
fi

# Determine max and min sizes and primer scheme for analysis
if [ $PRIMERS = "Sanger" ] ; then
    MAX=1250
    MIN=750
    PRIMERSCHEME=SARS-CoV-2/V1
elif [ $PRIMERS = "Artic" ] ; then
    MAX=700
    MIN=400
    PRIMERSCHEME=SARS-CoV-2/V3
else
    printf "${RED}ERROR:${NC} Unrecognised primer scheme\n\n"
    Help
    exit
fi

#Ensure the barcode list is unique
if [ $( awk -F"\t" '{print $1}' $SAMPLES | sort | uniq | wc -l ) != $( <$SAMPLES wc -l ) ] ; then
    printf "${RED}ERROR:${NC} Duplicate barcodes present in $SAMPLES file\n"
    exit
fi

#Ensure the Sample names are unique
if [ $( awk -F"\t" '{print $2}' $SAMPLES | sort | uniq | wc -l ) != $( cat $SAMPLES | wc -l ) ] ; then
    printf "${RED}ERROR:${NC} Duplicate samples present in $SAMPLES file\n"
    exit
fi

echo -e "${GREEN}CHECKED:${NC}All required programs, files and locations are present\n\n"

#==============THE SCRIPT==================Se   
#Run the guppy barcoder to demultiplex
if [ $GUPPYBARCODER = 1 ] ; then
    printf "${BLUE}1 - Running the guppy_barcoder to demultiplex the FASTQ files.${NC} \n\n"
    guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $FASTQ_OUT --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg
else
    printf "###### ${GREEN}Skipping guppy_barcoder step${NC} ######\n\n"
fi

#Move to the correct directory
cd $FASTQ_OUT

if [ $GUPPYPLEX = 1 ] ; then
    printf "${BLUE}2 - Combining demultiplexed files into a single fastq and excluding based on size.${NC}\n\n"
    #Run through the demuxed barcodes to combine into single fastq files
    for FILE in barcode[0-9]*; do
        artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $OUTPUT/fastq/$FILE --prefix $RUNNAME
    done
else
    printf "###### ${GREEN}Skipping guppyplex step${NC} ######\n\n"
fi


if [ $ARTICMINION = 1 ] ; then
    printf "${BLUE}3 - Importing samplenames and processing with artic minion command.${NC}\n\n"
    #Import samplenames / barcodes and then run the processing
    for FILE in $RUNNAME\_barcode[0-9]*.fastq ; do 
        BARCODE=`echo $FILE | sed s/$RUNNAME// | sed s/_barcode// | sed s/.fastq//`
        SAMPLENAME=`awk -F"\t" '$1=='$BARCODE' {print $2}' $SAMPLES | tr -d '\r' ` #tr to remove returns
        SAMPLENAME=${SAMPLENAME:- "Unknown_BC$BARCODE"}
        
        #Ensure fastq file is not empty
        FASTQLENGTH=`wc -l $FILE | sed s/$FILE//`  
        
        #Remove files that have not got enough data
        if [ $FASTQLENGTH -gt 200 ] ; then
            echo -e "\n\n Processing barcode number $BARCODE, Sample $SAMPLENAME \n"
            #Run processing scheme
            if [ $MEDAKA = 1 ] ; then
                printf "${GREEN}Using Medaka pipeline${NC}\n"
                artic minion --medaka --normalise 200 --threads 4 --scheme-directory $PRIMERDIR --read-file $FILE nCoV-2019/V1 $SAMPLENAME
            else
                printf "${GREEN}Using Nanopolish pipeline${NC}\n"
                artic minion --normalise 200 --threads 4 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $FAST5RAW --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME
            fi
            
            #Make samplename subdirectory if required
            if [ ! -d $SAMPLENAME ] ; then
                mkdir $SAMPLENAME
                echo "Making subdirectory $SAMPLENAME"
            fi
            
            #Move output from what was produced into its own directory
            find ./ -type f -name "$SAMPLENAME*" | xargs -I '{}'  mv {} "$SAMPLENAME"/
            #Move directory to another level for clarity
            mv $SAMPLENAME/ $FINAL_OUT/$SAMPLENAME
        else
            printf "${RED}ERROR:${NC} Too few reads in File $FILE (BC $BARCODE, Sample $SAMPLENAME).
            Aborting processing this file\n"
        fi
    done
else
    echo "Skipping artic minion step"
fi

echo -e "\nAll process completed successfully. \n\n"
exit

######################## NOTES ##########
#Once the above is complete

# Concatenate and tidy up the names of the fasta headers
cat *.fasta | sed 's/\/ARTIC\/nanopolish MN908947.3/\/2020/'| sed 's/>/>SARS-CoV-2\/human\/Zambia\//' > out.fasta

#Bug catching - list all defined variables
( set -o posix ; set ) | less

#Works:
artic minion --normalise 200 --threads 4 --scheme-directory /home/dan/Work/Computer/GitRepos/artic-ncov2019/primer-schemes/ --read-file C1/fastq/C1_barcode01.fastq --fast5-directory /media/dan/Master/WGS/C1_2020-09-09/ --sequencing-summary /media/dan/Master/WGS/C1_2020-09-09/sequencing_summary_FAO22552_3fe407a1.txt SARS-CoV-2/V3 29
