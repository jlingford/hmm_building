#!/usr/bin/bash

# make sure famsa binary is available on PATH
#
# famsa usage:
# famsa [options] infile outfile
#
# see `famsa -h`
#
# seems like each file needs to be generated with separate runs... won't output all in one go
#

# positional args
infaa=$1    # WARN: expects `.faa` extension
CPU=${2:-0} # default to 0 cores if not set

outaln=${infaa%.*}.famsa.aln
guidetree=${infaa%.*}.guidetree.dnd
log="$(dirname "$infaa")/famsa.log"

########### Just get the MSA ############

# 1. align sequences
famsa -v -t $CPU $infaa $outaln 2>$log

####### Saving tree and aln #############

# 1. export a single linkage guide tree to the Newick format
# NOTE: slowest step, but next steps should be faster
famsa -gt sl -gt_export -v -t $CPU $infaa $guidetree 2>$log

# 2. align sequences with the previously generated guide tree
# WARN: this step crashes?? report to github?
# famsa -gt import $guidetree -v -t $CPU $infaa ${outfile}.famsa.aln

######### OPTIONAL OUTPUTS ##############
## WARN: files can be very large
##
## optional: export a distance matrix to CSV format (square)
#famsa -dist_export -pid -square_matrix -v -t $CPU $infaa ${outfile}.dist.csv
#
## optional: export a pairwise identity (PID) matrix to the CSV format (square)
#famsa -dist_export -pid -square_matrix -v -t $CPU $infaa ${outfile}.pid.csv
