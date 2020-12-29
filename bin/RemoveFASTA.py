#!/usr/bin/python3.8
# ----------------------------------------------------------------------
# Script to selectively remove fasta headers and subsequent sequence from a file
# Version 0.1 29/12/2020
# Written by Daniel Bridges (danieljbridges@gmail.com)
#
#Changelog
# ----------------------------------------------------------------------

import sys

#Import the command line arguments into variables
fasta_path = str(sys.argv[1])
search = ">"+str(sys.argv[2])
fasta_output = str(sys.argv[1])+".final"

#Identify all lines to be retained into a list
final= []
match = 0

with open(fasta_path, "r") as fasta:
    for l in fasta:
        #Identify fasta header
        if l[0] == ">":
            #Identify fasta with the search string
            if l.startswith(search) :
                match = 1
            else :
            #Retain those without the string 
                match = 0
                final.append(l)
        elif match == 0 :
            final.append(l)

#Open outputfile and read in elements of the list above
outputfile = open(fasta_output, "w")
for line in final:
    outputfile.write(line)
outputfile.close()  
