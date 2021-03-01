#!/bin/bash
# ----------------------------------------------------------------------
# Dans script to synchronise between the MinIT and the laptop to allow RAMPART to run locally
# Version 1.0 09/09/20
# ----------------------------------------------------------------------
#Pass any CLI options
MINITPATH=$1
LOCALPATH=$2

MINITIP='10.42.0.1'
OPTIONS='-av -e'

function SYNCHRONISE ()
{
    #build rsync command
    rsync $OPTIONS ssh minit@$MINITIP:/$MINITPATH $LOCALPATH
}

function CHECK ()
{
if [[ -z $MINITPATH || -z $LOCALPATH ]]; then
  echo -e "\nYou have not supplied enough command line arguments (see usage above)."
  exit 1
fi
}
#==============THE SCRIPT==================

printf "This script will synchronise data from the MinIT (assumed to be at $MINITIP) to the local Laptop.
Usage:\n
MinIT2Laptop.sh <PathToMinIT> <LocalPATH>\n

e.g. 

MinIT2Laptop.sh data/C2_15-09-2020/C2_15-09-2020/20200916_1110_MN21231_FAO22552_20e8fe978fe97/ /home/dan/WGS/RawData/C2_15-09-2020/
 

Hit CTRL-C if this is not what you want to run....."

sleep 1

#Ensure command line options supplied
CHECK

#Synchronise the folders
SYNCHRONISE

exit
