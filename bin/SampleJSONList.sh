#!/bin/bash
# ----------------------------------------------------------------------
# Script to generate JSOn output of Barcode and samples
# Version 0.1 04/11/2020
# Written by Daniel Bridges (danieljbridges@gmail.com)
# ----------------------------------------------------------------------
set -e #exit whenever a command exits with a non zero status
set -u #treat undefined variables as errors
set -o pipefail #pipe will be considered successful if all the commands are executed without errors

#==============FUNCTIONS==================

function Help {
   # Display Help
   echo -e "This script takes a tab delimited file listing the barcode (2 digits) and the sample name and generates a JSON file called run_configuration.json for use by rampart.
   Syntax: $(basename "$0") [-h|i|o|r]
      
    Required arguments:
   -i      Tab delimited file containing the assigned barcode (2 digit) and sample name
   -o      Directory location where the output file should be stored
   -r      Runname for rampart heading
   
    Optional arguments:
   -h      Print this help"
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
while getopts 'hi:o:r:' opt; do
  case "$opt" in
    i)
      INPUTFILE=$OPTARG
      ;;
    h)
      Help
      exit
      ;;
    o)
      OUTPUTDIR=$OPTARG
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



#==============THE SCRIPT==================
#Give default entries
#INPUTFILE=${INPUTFILE:-0)
#OUTPUTDIR=${OUTPUTDIR:-0)
#RUNNAME=${RUNANME:-0)

#IF [ $INPUTFILE =0 |]


#Ensure files / directories are present
present $INPUTFILE "f"
present $OUTPUTDIR "d"

OUTPUTFILE="$OUTPUTDIR/run_configuration.json"

#Start the JSON File
printf "{
  \"title\": \"Run $RUNNAME\",
  \"basecalledPath\": \"fastq_pass\",
  \"samples\": [" > $OUTPUTFILE

#Loop through each barcode to create the JSON file
NUMBARCODES=`cat $INPUTFILE | wc -l`

COUNT=1
until [ $NUMBARCODES -lt $COUNT ]; do
    #Pull out the data sequentially by row (FNR)
    SAMPLE=`awk -F"\t" 'FNR == '$COUNT' {print $2}' $INPUTFILE | tr -d '\r'`
    BARCODE=`awk -F"\t" 'FNR == '$COUNT' {print $1}' $INPUTFILE`
    printf "    {
      \"name\": \"$SAMPLE ($BARCODE)\",
      \"description\": \"\",
      \"barcodes\": [ \"NB$BARCODE\" ]" >> $OUTPUTFILE
    if [ $NUMBARCODES == $COUNT ]; then
        echo -e "    }" >> $OUTPUTFILE #If last entry then there shouldn't be a comma
    else
        echo -e "    }," >> $OUTPUTFILE 
    fi
    #Increment the count
    let COUNT=COUNT+1
    printf "Sample = $SAMPLE, Barcode = $BARCODE\n"
done

#Finish off the JSON File
printf "  ]
}" >> $OUTPUTFILE

