#!/bin/bash
# ----------------------------------------------------------------------
# Script to generate JSOn output of Barcode and samples
# Version 0.1 04/11/2020
# Written by Daniel Bridges (danieljbridges@gmail.com)
# ----------------------------------------------------------------------
set -e #exit whenever a command exits with a non zero status
set -u #treat undefined variables as errors
set -o pipefail #pipe will be considered successful if all the commands are executed without errors


#Recognise any CLI options
while getopts ':i:o:r:' opt; do
  case "$opt" in
    i)
      INPUTFILE=$OPTARG
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

#==============FUNCTIONS==================

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

#==============THE SCRIPT==================

present $INPUTFILE "f"
present $OUTPUTDIR "d"

OUTPUTFILE="$OUTPUTDIR/run_configuration.json"

#Start the JSON File
echo -e "{
  \"title\": \"Run $RUNNAME\",
  \"basecalledPath\": \"fastq_pass\",
  \"samples\": [" > $OUTPUTFILE

#Loop through each barcode to create the JSON file
NUMBARCODES=`cat $INPUTFILE | wc -l`

COUNT=1
until [ $NUMBARCODES -lt $COUNT ]; do
    #Pull out the data sequentially by row (FNR)
    SAMPLE=`awk -F"\t" 'FNR == '$COUNT' {print $2}' $INPUTFILE`
    BARCODE=`awk -F"\t" 'FNR == '$COUNT' {print $1}' $INPUTFILE`
    echo -e "    {
      \"name\": \"$SAMPLE\",
      \"description\": \"\",
      \"barcodes\": [ \"NB$BARCODE\" ]" >> $OUTPUTFILE
    if [ $NUMBARCODES == $COUNT ]; then
        echo -e "    }" >> $OUTPUTFILE #If last entry then there shouldn't be a comma
    else
        echo -e "    }," >> $OUTPUTFILE 
    fi
    #Increment the count
    let COUNT=COUNT+1
done

#Finish off the JSON File
echo -e "  ]
}" >> $OUTPUTFILE

