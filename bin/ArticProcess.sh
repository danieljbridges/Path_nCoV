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

VERSION="0.7"
VERDATE="2021-06-06"
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
   
   Syntax: $(basename "$0") [-123456hv|b|r|m|p|s]
   Version: $VERSION, $VERDATE
   
    Mandatory args:
    -b      Basefolder for all inputs and outputs - please use an absolute path
    -r      Runname to prefix output files, folder names etc

    Optional args:
    -h      Print this help
    -m      Process using the medaka pipeline rather than with nanopolish
    -p      Primers used during PCR [Sanger | Artic | Midnight]
    -s      Location directory of the primer scheme used
    -v      Version

    Optional steps to omit:
    -1      Omit Step 1 (guppy_barcoder)
    -2      Omit Step 2 (guppyplex)
    -3      Omit Step 3 (artic minion)
    -4      Omit Step 4 (Consensus sequences)
    -5      Omit Step 5 (Sequencing statistics)
    -6      Omit Step 6 (Generate JSON file for RAMPART)
    -7      Omit Step 7 (All sequences and filtered list of sequences)
    -8      Omit Step 8 (ncov-tools QC pipeline)
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
        printf "${RED}ERROR:${NC} $2 not present. Please install the following and try again:\n"
        printf "   $3\n\n"
        printf "${RED}Exiting script${NC}\n\n"
        exit
    fi
}

function GETCOLMID {
awk -F"$1" 'NR==1{ for (i=1;i<=NF;i++) if ($i == "'$2'") print i }' $3
}

function change_conda {
    CONDA=`which conda | sed s'/\/condabin//' | sed s'/bin\///' | sed s'/\/conda//'`
    source $CONDA/etc/profile.d/conda.sh 
    set +eu
    conda activate $1
}

#==============Recognise any CLI options==================
while getopts 'vhp:b:r:s:m12345678' opt; do
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
    7)
      S7=0
      ;;
    8)
      S8=0
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
printf "\n${ORANGE}##### Running `basename $0` Version $VERSION ##### ${NC}\n"
printf "${GREEN}Checking all files, folders and data are present${NC}\n"
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
S8=${S8:- 1}
MEDAKA=${MEDAKA:- 0}            # 3
PRIMERDIR=${PRIMERDIR:- 1}      # 3
PRIMERS=${PRIMERS:- 1}          # 2, 3
ARTIC_OUT="$BASEFOLDER/3_Artic_Output/$RUNNAME"
RAWDATADIR="$BASEFOLDER/1_Raw/$RUNNAME"

#################################################
#Make directories for the demultiplexed, combined and size selected fast_q files
check_mkdir "$ARTIC_OUT/fastq"
check_mkdir "$ARTIC_OUT/processed"
present $RAWDATADIR "d"

ALLSEQ="allsequences.fasta" #Name of file for all sequences to output to

#Determine if the necessary files and directories and arguments etc exist for each step
printf ""
if [ $S1 = 1 ] ; then
    check_package "guppy_barcoder" "guppy" "guppy (from ONT website)\n"
    FASTQRAW="$RAWDATADIR/fastq_pass"
    present $FASTQRAW "d"
fi

if [ $S2 = 1 ] || [ $S3 = 1 ]; then
    check_package "artic" "artic" "artic"
fi

if [ $S2 = 1 ] || [ $S3 = 1 ] || [ $S8 = 1 ]; then
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
    elif [ $PRIMERS = "Midnight" ] ; then
        MAX=950
        MIN=1450
        PRIMERSCHEME=SARS-CoV-2/V2
    else
        Help
        printf "${RED}ERROR:${NC} Unrecognised primer scheme\n\n"
        exit
    fi
fi

if [ $S3 = 1 ] || [ $S8 = 1 ]; then
    #Check primer scheme supplied
    check_var $PRIMERDIR "-s (primer scheme location)"
    PRIMERDIR=`realpath $PRIMERDIR`
fi


if [ $S3 = 1 ] ; then
    SEQUENCINGSUMMARY=$(find $RAWDATADIR -name 'sequencing_summary*' -not -name '*.tmp')
    check_var ${SEQUENCINGSUMMARY:- 1} "Sequencing Summary"
fi

if [ $S3 = 1 ] || [ $S6 = 1 ] || [ $S8 = 1 ] ; then
    #Identify csv file
    SAMPLEFILE="$BASEFOLDER/2_SampleList_and_Rampart/Samples_Sequenced.csv"
    present $SAMPLEFILE "f"
fi


if [ $S3 = 1 ] || [ $S6 = 1 ] ; then
    #Check if all of the seqIDs are unique in the samples file
    readarray -t IDS < <(awk -F"," 'NR>1 {if ($8) print$8}' $SAMPLEFILE | sort | uniq -d)
    if [ ${#IDS[@]} -gt 0 ]; then
        printf "${RED}ERROR:${NC} The following duplicate seqIDs exist in $SAMPLEFILE. 
        Please correct so that every seqID (column 8) is unique to the entire list (empty records are not included).\n"
        printf '%s\n' "${IDS[@]}"
        exit
    fi
    
    #Read in Samples
    readarray -t SAMPLES < <(awk -F"," '$9 ~ /'$RUNNAME'/ {print$8}' $SAMPLEFILE)
    if (( ${#SAMPLES[@]} == 0 )); then
        printf "\n${RED}ERROR:${NC} SAMPLES array is empty \nExiting script\n\n"
        exit
    fi
    
    #Read in barcodes
    readarray -t BARCODES < <(awk -F"," '$9 ~ /'$RUNNAME'/ {print$7}' $SAMPLEFILE)
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
    printf "${GREEN}CHECKED:${NC}No errors in sample list.\n"
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
    printf "guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files barcode_arrs_nb96.cfg${NC}\n" | tee -a "${RUNLOG}1.log"
    
    guppy_barcoder --require_barcodes_both_ends -i $FASTQRAW -s $ARTIC_OUT/fastq --arrangements_files "barcode_arrs_nb96.cfg" | tee -a "${RUNLOG}1.log"
    printf "\n###### ${GREEN}Step 1: guppy_barcoder completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 1: Skipping guppy_barcoder step${NC} ######\n\n"
fi

#STEP 2: Combine all identical barcode reads into same fasta and exclude sizes
if [ $S2 = 1 ] ; then
    printf "\n###### ${BLUE}Step 2: Combining demultiplexed files into a single fastq and excluding based on size.${NC} ######\n\n" | tee "${RUNLOG}2.log"
    readarray -t S2DIRS < <(find $ARTIC_OUT/fastq -type d -name 'barcode[0-9]*')
    
    #Change directory as unable to redirect output from guppyplex
    cd $ARTIC_OUT/fastq
    
    #Run through the demuxed barcodes to combine into single fastq files
    for S2DIR in "${S2DIRS[@]}" ; do
        printf "${LG}artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME ${NC}\n" | tee -a "${RUNLOG}2.log"
       artic guppyplex --skip-quality-check --min-length $MIN --max-length $MAX --directory $S2DIR --prefix $RUNNAME | tee -a "${RUNLOG}2.log"
    done
    echo -e "\n###### ${GREEN}Step 2: guppyplex completed. ${NC} ######\n\n"
else
    printf "###### ${GREEN}Step 2: Skipping guppyplex step${NC} ######\n\n"
fi

#STEP 3: Processing with artic minion pipeline
if [ $S3 = 1 ] ; then
    printf "\n###### ${BLUE}Step 3: Importing samplenames and processing with artic minion command.${NC} ######\n\n" | tee "${RUNLOG}3.log"
    
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
                    printf "\n${GREEN}Using Medaka pipeline${NC}\n"
                    printf "${LG}artic minion --medaka --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE $PRIMERSCHEME $SAMPLENAME ${NC}\n" | tee -a "${RUNLOG}3.log"
                    artic minion --medaka --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${RUNLOG}3.log"
                    # 2>&1 redirects stderr to stdout
                else
                    printf "\n${GREEN}Using Nanopolish pipeline${NC}\n"
                    printf "${LG}artic minion --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME${NC}\n" | tee -a "${RUNLOG}3.log"
                    artic minion --normalise 200 --threads 24 --scheme-directory $PRIMERDIR --read-file $FILE --fast5-directory $RAWDATADIR --sequencing-summary $SEQUENCINGSUMMARY $PRIMERSCHEME $SAMPLENAME 2>&1 | tee -a "${RUNLOG}3.log"
                fi
                
                #Make samplename subdirectory if required
                check_mkdir $SAMPLENAME
                #Move output from what was produced into its own directory
                find ./ -type f -name "$SAMPLENAME*" | xargs -I '{}'  mv {} "$SAMPLENAME"/
                #Move directory to another level for clarity
                mv $ARTIC_OUT/fastq/$SAMPLENAME $ARTIC_OUT/processed/$SAMPLENAME
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
    printf "\n###### ${BLUE}Step 4: Concatenating consensus sequences from artic pipeline.${NC} ######\n\n"
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
    
    echo -e "\n###### ${GREEN}Step 4: Consensus sequences compiled. ${NC} ######\n\n"
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
    change_conda base
    check_package "nextclade" "nextclade package" "nextclade (https://www.npmjs.com/package/@neherlab/nextclade)"
    #run nextclade
    if [ $PRIMERDIR == 1 ] ; then
        nextclade -i $ALLSEQ -c nextclade.csv -o nextclade.json cd2>&1 | tee "${LOGFOLDER}nextclade.log"
    else
        PRIMERFILE="$PRIMERDIR/SARS-CoV-2/nextclade_primers.csv"
        printf "Using custom primer file ($PRIMERFILE)\n"
        present $PRIMERFILE "f"
        nextclade -i $ALLSEQ -c nextclade.csv -o nextclade.json --input-pcr-primers $PRIMERFILE cd2>&1 | tee "${LOGFOLDER}nextclade.log"
    fi
    mv nextclade* intermediates/
    
    printf "\n###### ${GREEN} Determining PANGO lineages ${NC} ######\n\n"
    change_conda pangolin
    check_package "pangolin" "pangolin environment" "pangolin"
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
    printf "\n###### ${BLUE}Step 7: Generate filtered fasta file of submittable entries. ${NC} ######\n\n"
    
    cd "$BASEFOLDER/5_GISAID"
    ALLSEQDATAFN="allsequencedata.csv"
    present $ALLSEQDATAFN "f"
    
    printf "\n###### ${GREEN} Filtering to remove sequences that are not submittable ${NC} ######\n\n"
    
    #Pull out the ID for specific column headers
    SUBCOL="$(GETCOLMID "\t" "Submittable" "$ALLSEQDATAFN")"
    SEQIDCOL="$(GETCOLMID "\t" "SeqID" "$ALLSEQDATAFN")"
        
    #Generate list of sequences
    awk -F"\t" -v i=$SUBCOL -v j=$SEQIDCOL '$i ~ /True/ {print$j}' allsequencedata.csv > samples.csv
    
    #Look for duplicates in the list (second check)
    DUPES=`cat samples.csv | sort | uniq -d`
    if [ -n "$DUPES" ] ; then
        printf "${RED} ERROR: Duplicate entries found in the list of SeqIDs to retain. \n Exiting script\n ${NC}\n"
        printf "$DUPES\n"
        exit
    fi
    
    #Filter the list of all sequences to only retain the correct ones
    seqkit grep -n -f samples.csv $ALLSEQ -o allsequences.submittable.fasta
    rm samples.csv
    
    MINCOV=90
    MINBREADTH=90
    
    printf "\n###### ${GREEN} Filtering to remove sequences with that are >${MINCOV} percent coverage and >${MINBREADTH}x average depth ${NC} ######\n\n"
    
    #Pull out all of the submittable SeqID entries that should be retained by searching column headers
    COVCOL="$(GETCOLMID "\t" "sequencing_depth_avg" "$ALLSEQDATAFN")"
    BREADTHCOL="$(GETCOLMID "\t" "coverage_breadth" "$ALLSEQDATAFN")"
    
    awk -F"\t" -v i=$COVCOL -v j=$BREADTHCOL -v k=$SUBCOL -v l=$SEQIDCOL '$i > 90 && $j > 90 && $k ~ /True/ {print$l}' allsequencedata.csv > samples.csv

    #Look for duplicates in the list (second check)
    DUPES=`cat samples.csv | sort | uniq -d`
    if [ -n "$DUPES" ] ; then
        printf "${RED} ERROR: Duplicate entries found in the list of SeqIDs to retain. \n Exiting script\n ${NC}\n"
        printf "$DUPES\n"
        exit
    fi
    
    #Filter the list of all sequences to only retain the correct ones
    seqkit grep -n -f samples.csv $ALLSEQ -o allsequences.submittable.qc.fasta
    rm samples.csv
else
    printf "###### ${GREEN}Step 7: Skipping generation of filtered fasta file${NC} ######\n\n"

fi

if [ $S8 = 1 ] ; then
    printf "###### ${BLUE}Step 8: Running ncov-tools QC pipeline${NC} ######\n\n"
    change_conda "ncov-qc"
    check_package "snakemake" "ncov-tools environment" "ncov-qc"
    QCDIR="$BASEFOLDER/6_QCAnalysis/$RUNNAME"
    check_mkdir $QCDIR
       
    cd "$BASEFOLDER/6_QCAnalysis/$RUNNAME"
    
    printf "Removing prior QC run details:\n"
    rm -rvf $QCDIR/*
    
    printf "\nGenerating config.yaml file\n"
    SEQIDCOL="$(GETCOLMID "," "SeqID" "$SAMPLEFILE")"
    CTCOL="$(GETCOLMID "," "Ctvalue" "$SAMPLEFILE")"
    SEQRUNCOL="$(GETCOLMID "," "SeqRun" "$SAMPLEFILE")"
    METADATAFN=Metadata.tsv
    awk -F"," -v i=$SEQRUNCOL -v k=$SEQIDCOL -v l=$CTCOL 'NR==1 {print "sample","ct"} $i ~ /'$RUNNAME'/ {print$k,$l}' OFS='\t' $SAMPLEFILE > $METADATAFN
    NTCLIST1=`awk -F"\t" '{print$1}' Metadata.tsv | grep NTC | sed -e 's/^/\\"/' | sed -e 's/$/\\"/' | awk -vORS=, '{print$1}'`
    NTCLIST=() #declare as an empty variable
    #ncov-tools can't deal with NTC that has no related bam file therefore if no bam remove from the list
    for NTC in "${NTCLIST1[@]}" ; do
        N=`echo $NTC | sed 's/"//g' | sed 's/,//g'`
        printf "N = $N\n"
        NTCDIR="${BASEFOLDER}/3_Artic_Output/$RUNNAME/processed/$N"
        printf "NTCDIR = $NTCDIR\n"
        if [ ! -d $NTCDIR ] ; then
            printf "Missing NTC directory for $N. Dropping from analysis\n"
        else
            NTCLIST+=($NTC)
        fi
    done
    
    YAMLFN="config.yaml"
    printf "Generating Metadata.tsv file\n"
    printf "data_root: $BASEFOLDER/3_Artic_Output/$RUNNAME/processed
run_name: $RUNNAME
reference_genome: ${PRIMERDIR}/${PRIMERSCHEME}/SARS-CoV-2.reference.fasta
platform: oxford-nanopore
primer_bed: ${PRIMERDIR}/${PRIMERSCHEME}/SARS-CoV-2.bed
bam_pattern: \"{data_root}/{sample}/{sample}.sorted.bam\"
consensus_pattern: \"{data_root}/{sample}/{sample}.consensus.fasta\"
variants_pattern: \"{data_root}/{sample}/{sample}.pass.vcf.gz\"
skip_empty_negatives: true
metadata: \"$METADATAFN\"
negative_control_samples: [ $NTCLIST ]\n" > $YAMLFN
    
    printf "Running snakemake report\n"
    snakemake --cores all -s ~/git/ncov-tools/workflow/Snakefile all_final_report
else
    printf "###### ${GREEN}Step 8: Skipping generation of QC stats${NC} ######\n\n"
fi
exit
