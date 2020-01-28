#!/usr/bin/env bash

set -u
set -o pipefail

################################################################################
# PRELIMINARIES

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=4
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
ConfigFile="$1"
bam="$2"
TheRef="$3"
SID="$4"

# Source required code & check files exist
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ConfigFile" "$bam"
CheckConfig "$ConfigFile" false false true || \
{ echo "Problem with $ConfigFile. Quitting." >&2 ; exit 1 ; }


# Some files we'll create
consensus="$SID"'_MinCov_'"$MinCov1"'_'"$MinCov2.fasta"
BaseFreqs="$SID$BaseFreqsSuffix"
InsertSizeCounts="$SID$InsertSizeCountsSuffix"
################################################################################

cp "$TheRef" .
TheRef=$(basename "$TheRef")

# Index
"$samtools" faidx "$TheRef" || \
{ echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
exit 1 ; }
"$samtools" index "$bam" 

# TODO: separate the mapping function out into mapping, and the bit below, then
# replace the bit below by its function.

# Sort the bam file. Thanks Nick Croucher!
"$samtools" view -b $samtoolsReadFlags -t "$TheRef".fai -o \
"$MapOutConversion1" "$bam" &&
"$samtools" sort -n "$MapOutConversion1" "$MapOutConversion2" &&
"$samtools" fixmate "$MapOutConversion2".bam "$MapOutConversion3" &&
"$samtools" sort "$MapOutConversion3" "$SID" &&
"$samtools" index "$SID.bam" || \
{ echo 'Failed to sort the bam file. Quitting.' >&2 ; exit 1 ; }

# Calculate the normalised insert size distribution.
#TODO: uncomment
#"$samtools" view "$bam" | awk '{if ($9 > 0) print $9}' > "$InsertSizes1"
#sort -n "$InsertSizes1" | uniq -c > "$InsertSizes2"
#InsertCount=$(awk '{sum+=$1} END {print sum}' "$InsertSizes1")
#awk '{print $2 "," $1 "," $1/'$InsertCount'}' "$InsertSizes2" > \
#"$InsertSizeCounts"

# Generate pileup
"$samtools" mpileup $mpileupOptions -f "$TheRef" "$SID.bam" > "$PileupFile" || \
{ echo 'Failed to generate pileup. Quitting.' >&2 ; exit 1 ; }

# Generate base frequencies and consensuses
"$Code_AnalysePileup" "$PileupFile" "$TheRef" > "$BaseFreqs" && \
"$Code_CallConsensus" "$BaseFreqs" "$MinCov1" "$MinCov2" > \
"$consensus" || \
{ echo 'Problem analysing the pileup or calling the consensus.' >&2 ; exit 1 ; }


