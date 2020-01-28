#!/usr/bin/env python2
from __future__ import print_function
#
## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
#
## Overview: takes a fastq file sorted by read name, e.g. with
## $ cat MyReads.fq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > out.fq
## and a file of read names sorted in the same way (one read name per line), and
## finds those reads, printing them to the screen suitable for redirection into
## another fastq file.
## Usage:
## $ FindNamedReadsInSortedFastq.py AllReadsSorted.fastq ReadNamesIwant.txt

# Import some modules we'll need.
import os.path
import sys
from optparse import OptionParser

# Define the arguments and options
parser = OptionParser()
parser.add_option("-v", "--invert", action="store_true", dest="invert",
default=False,
help="find all reads EXCEPT those named.")
(options, args) = parser.parse_args()
invert = options.invert

# Check that this script was called from the command line with two arguments.
if len(args) != 2:
  sys.stderr.write('Incorrect number of arguments given. Correct usage:\n '+\
  sys.argv[0] +' SortedFastqFile SortedFileOfReadNames\nQuitting\n')
  exit(1)
FastqFile = args[0]
FileOfReadNames = args[1]

# Check that all arguments exist and are files
for DataFile in [FastqFile,FileOfReadNames]:
  if not os.path.isfile(DataFile):
    sys.stderr.write(DataFile +' does not exist or is not a file. Quitting.\n')
    exit(1)

# Read in the read names
ReadNames = []
with open(FileOfReadNames, 'r') as f:
  for line in f:
    ReadNames.append(line.strip())
NumReadsToFind = len(ReadNames)

# Read through the fastq file...
NumNamedReadsFound = 0
ThisIsAReadWeWant = False
with open(FastqFile, 'r') as f:
  for LineNumberMin1, line in enumerate(f):

    # Those lines that are not read names: print if appropriate, and continue.
    if LineNumberMin1 % 4 != 0:
      if ThisIsAReadWeWant:
        print(line.rstrip())
      continue

    # Henceforth we're on a line that's a read name (lines 1, 5, 9, 13...)

    # Check if we've found all the named reads already.
    if NumNamedReadsFound == NumReadsToFind:
      if invert:
        ThisIsAReadWeWant = True
        print(line.rstrip())
        continue
      else:
        break

    # Check it begins with an @ symbol
    if not (len(line) > 0 and line[0] == '@'):
      sys.stderr.write('Unexpected fastq format for ' +FastqFile+ ': lines' +\
      ' 1, 5, 9, 13... are expected to start with an @ symbol. Quitting.\n')
      exit(1)

    # See if this is a read we want: a named read if invert = False,
    # or an unnamed read if invert = True.
    if line[1:].rstrip() == ReadNames[NumNamedReadsFound]:
      NumNamedReadsFound += 1
      ThisIsAReadWeWant = not invert
    else:
      ThisIsAReadWeWant = invert
    if ThisIsAReadWeWant:
      print(line.rstrip())

# Check all reads were found
if NumNamedReadsFound < NumReadsToFind:
  sys.stderr.write(ReadNames[NumNamedReadsFound] +' was not found in ' +\
  FastqFile +'. Either it is truly missing, or ' + FastqFile +\
  ' is not sorted, or ' +FileOfReadNames +' is not sorted, or'+FileOfReadNames+\
  ' contains a multiply specified read.\nQuitting.\n')
  exit(1)

