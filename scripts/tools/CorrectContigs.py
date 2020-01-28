#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk and Francois Blanquart who
## wrote the original version of this script in R. Thanks to Tanya Golubchick
## for suggesting discarding contig sequence not contained in blast hits.
## Acknowledgement: this was written thanks to ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Checks whether contigs appear to need correction by
analysing a comma-separated blast file with variables
qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send
resulting from blasting contigs to a set of reference sequences. By default, we
exit with an error if correction is needed - if contigs have multiple hits
(ignoring any hit fully inside another) or hits are in the reverse direction.
Optionally, the contigs themselves can be supplied and correction is attempted:
cutting the contigs where they have multiple hits and/or reverse-complementing
them where hits are in the reverse direction.
'''

import argparse
import os
import sys
import copy
from Bio import SeqIO
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BlastFile', type=File)
parser.add_argument('-C', '--contigs', type=File, \
help='The fasta file of contigs.')
parser.add_argument('-O', '--out-file', help='''The file to which corrected
contigs will be written. Specifying "-" will print output to stdout. If
correction is not required, this file will not be created.''')
parser.add_argument('--overwrite', action='store_true',
help='If the out file exists already, overwrite it instead of stopping.')
parser.add_argument('-F', '--min-hit-frac', type=float, help='''Only relevant if
we're just checking whether correction is needed, not if we're actually doing
the correction. Use this option to specify a minimum fraction of a contig's
length that its hit must cover: below this we assume correction is needed
(either discarding the whole contig or removing the non-blasting segment). The
default is 0.9.''', default=0.9)
parser.add_argument('-K', '--keep-non-hits', action='store_true', help='''Only
relevant if contig correction is being performed. By default we discard non-hit
regions: stretches of contig sequence not contained inside any blast hit. With
this option, such regions are kept. If there is such a region in between two
hits, it will be duplicated so that it appears as part of the new contig we
create for the hit to the left, and also for the one to the right. e.g. say your
contig looked like this: 1111ACGT2222 where 1111 denotes sequence covered by one
blast hit, 2222 is sequence covered by a second blast hit, and the ACGT in
between did not blast. By default we would produce the two contigs 1111 and
2222; with this option we produce instead 1111ACGT and ACGT2222. See also
--dont-duplicate which interacts with this option.''')
parser.add_argument('-D', '--dont-duplicate', action='store_true', \
help='''Only relevant if contig correction is being performed. When two hits
partially overlap, we cut to create two new contigs; by default the overlap
(i.e. the sequence contained in both hits) is duplicated to be present in both
of the new contigs. With this option, we instead cut half way through the
overlap, dividing the overlapping region between the two new contigs created by
the cut. e.g. say your contig looked like this: 1111ACGT2222 where 1111 denotes
sequence uniquely contained in one blast hit, 2222 is sequence uniquely
contained in a second blast hit, and the ACGT in between is contained in both
blast hits. By default we would produce the two contigs 1111ACGT and ACGT2222;
with this option we produce instead 1111AC and GT2222. The same results would be
given if the ACGT was in neither hit (instead of both) and you are using
--keep-non-hits: 1111ACGT and ACGT2222 without this option, 1111AC and GT2222
with it.''')
args = parser.parse_args()

# The --contigs and --out-file flags need each other.
MakeCorrections = args.contigs != None
if MakeCorrections and args.out_file == None:
  print('The --contigs and --out-file flags need each other: use both or', \
  'neither. Quitting.', file=sys.stderr)
  exit(1)

# If we're writing output: use stdout if desired, otherwise check the output
# file doesn't exist and is writable.
if MakeCorrections:
  if args.out_file == "-":
    args.out_file = sys.stdout
  elif not args.overwrite and os.path.isfile(args.out_file):
    print(args.out_file, 'exists already. Move, rename or delete it, and try',
    'again. Quitting.', file=sys.stderr)
    exit(1)

# Check the min hit frac is in (0,1)
if not (0 < args.min_hit_frac < 1):
  print('The --min-hit-frac value should be greater than 0 and less than 1.', \
  'Quitting.', file=sys.stderr)
  exit(1)

# columns of the blast output are:
# qseqid means Query Seq-id
# sseqid means Subject Seq-id
# evalue means Expect value
# pident means Percentage of identical matches
# qlen means Query sequence length
# qstart means Start of alignment in query
# qend means End of alignment in query
# sstart means Start of alignment in subject
# send means End of alignment in subject

# Read in the blast hits.
HitDict = collections.OrderedDict()
with open(args.BlastFile) as f:
  for line in f:
    try:
      qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send = \
      line.split(',')
    except ValueError:
      print('The following line in ', args.BlastFile, ' does not contain 9 ', \
      'columns:\n', line, 'Quitting.', sep='', file=sys.stderr)
      exit(1)
    try:
      evalue, pident, qlen, qstart, qend, sstart, send = float(evalue), \
      float(pident), int(qlen), int(qstart), int(qend), int(sstart), int(send)
    except ValueError:
      print('Could not understand columns 3-9 (evalue, pident, qlen, qstart,' ,\
      'qend, sstart, send) on line\n', line, 'in ', args.BlastFile, ' as',\
      'floats & ints. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if qstart >= qend:
      print('qstart greater than or equal to qend on line\n', line, 'in ', \
      args.BlastFile + '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if sstart == send:
      print('sstart = send on line\n', line, 'in ', args.BlastFile + \
      '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if min(qstart, sstart, send) < 1:
      print('Non-positive value of qstart, sstart or send on line\n', line, \
      'in ', args.BlastFile + '. Unexpected. Quitting.', sep='', \
      file=sys.stderr)
      exit(1)
    if qend > qlen:
      print('qend greater than qlen on line\n', line, 'in ', args.BlastFile + \
      '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    hit = [qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send]
    if qseqid in HitDict:
      HitDict[qseqid].append(hit)
    else:
      HitDict[qseqid] = [hit]

# Quit if no hits.
if len(HitDict) == 0:
  print(args.BlastFile, 'contains no hits. Quitting.', file=sys.stderr)
  exit(1)

CorrectionsNeeded = False
for contig, hits in HitDict.items():

  # Where a contig has one hit contained entirely inside another, remove the
  # sub-hit.
  if len(hits) > 1:
    SubHitIndices = []
    for i, hit1 in enumerate(hits):
      if i in SubHitIndices:
        continue
      qstart1, qend1 = hit1[5:7]
      qmin1, qmax1 = min(qstart1, qend1), max(qstart1, qend1)
      for j, hit2 in enumerate(hits[i+1:]):
        j = i + j + 1 # offset & zero-based indexing
        if j in SubHitIndices or i in SubHitIndices:
          continue
        qstart2, qend2 = hit2[5:7]
        qmin2, qmax2 = min(qstart2, qend2), max(qstart2, qend2)
        # If the hits start and end at the same place, only one should go:
        if qmin1 == qmin2 and qmax1 == qmax2:
          SubHitIndices.append(j)
        elif qmin1 <= qmin2 and qmax1 >= qmax2:
          SubHitIndices.append(j)
        elif qmin1 >= qmin2 and qmax1 <= qmax2:
          SubHitIndices.append(i)
    assert len(SubHitIndices) == len(set(SubHitIndices)), \
    'Internal error removing sub-hits. Please report to Chris Wymant.'
    for i in sorted(SubHitIndices, reverse=True):
      del HitDict[contig][i]
    hits = HitDict[contig]

    # If we're only checking (not correcting) and there are multiple hits or
    # a reverse hit or a too small hit, quit.
    if len(hits) > 1:
      if not MakeCorrections:
        print('Contig correction required (contig', contig, 'has multiple hits,',\
        'after removing any that are fully contained within another). Quitting.',\
        file=sys.stderr)
        exit(1)
      CorrectionsNeeded = True

  FirstHit = hits[0]
  qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send = FirstHit

  # If there's a bit of the contig that doesn't blast and the user wants to
  # discard such bits, correction is needed. However if we're just checking 
  # whether correction is needed (i.e. not actually correcting the contigs), 
  # don't quit just yet: the --min-hit-frac option is designed to allow some 
  # flexibility so that we don't report that correction is necessary because 
  # 0.001% of the contig fails blasting.
  HitLength = qend - qstart + 1
  if HitLength < qlen and not args.keep_non_hits and MakeCorrections:
    CorrectionsNeeded = True

  # Quit if the hit fraction is too small, if desired.
  HitFrac = float(HitLength) / qlen
  if not MakeCorrections and HitFrac < args.min_hit_frac:
    print('Contig correction required (the hit\n', \
    ' '.join(map(str, FirstHit)), '\nfor contig', contig, 'has a hit fraction',\
    HitFrac, "which is below the --min-hit-frac value of", \
    str(args.min_hit_frac) + "). Quitting.", file=sys.stderr)
    exit(1)

  # If reverse complementation is needed, correction is needed. Quit if desired.
  if sstart > send:
    if not MakeCorrections:
      print('Contig correction required (the hit\n', \
      ' '.join(map(str, FirstHit)), '\nfor contig', contig + \
      "is in the reverse direction). Quitting.", file=sys.stderr)
      exit(1)
    CorrectionsNeeded = True

# If no corrections are needed, quit successfully.
if not CorrectionsNeeded:
  print('No contig correction needed for', args.BlastFile + '.')
  exit(0)

# Read in the contigs
ContigDict = collections.OrderedDict()
for seq in SeqIO.parse(open(args.contigs),'fasta'):
  if seq.id in ContigDict:
    print('Encountered sequence', seq.id, 'a second time in', args.contigs+\
    'Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)
  ContigDict[seq.id] = seq

# Check we have a sequence for each hit
UnknownHits = [hit for hit in HitDict.keys() if not hit in ContigDict.keys()]
if len(UnknownHits) != 0:
  print('The following hits in', args.BlastFile, 'do not have a corresponding',\
  'sequence in', args.contigs +':\n', ' '.join(UnknownHits) + \
  '\nQuitting.', file=sys.stderr)
  exit(1)

OutSeqs = []
for ContigName, hits in HitDict.items():

  seq = ContigDict[ContigName]
  SeqLength = len(seq.seq)

  # Check all hits for a contig report the same length that we observe.
  for hit in hits:
    qlen = hit[4]
    if qlen != SeqLength:
      print(args.BlastFile, 'has qlen =', qlen, 'for contig', ContigName, \
      'but that contig has length', SeqLength, 'in', args.contigs + \
      '. Quitting.', file=sys.stderr)
      exit(1)

  # If a contig has only one hit, easy: trim off the bits not covered by the
  # hit if desired, and reverse complement if needed.
  NumHits = len(hits)
  if NumHits == 1:
    if not args.keep_non_hits:
      qstart, qend = hits[0][5:7]
      seq.seq = seq.seq[qstart-1:qend]
    sstart, send = hits[0][7:9]
    if sstart > send:
      seq.seq = seq.seq.reverse_complement()
    seq.id += '_BlastsTo_' + str(min(sstart, send)) + '-' + \
    str(max(sstart, send))
    OutSeqs.append(seq)

  else:

    # Sort hits by their start point.
    StartPoints = [hit[5] for hit in hits]
    assert len(StartPoints) == len(set(StartPoints)), \
    'Internal error cutting contigs. Please report to Chris Wymant.'
    hits = sorted(hits, key=lambda x:x[5])

    for i, hit in enumerate(hits):

      ThisStart, ThisEnd = hit[5:7]
      if not args.dont_duplicate:

        # For this block, we want to keep sequence that's not inside a hit, and
        # duplicate sequence that's between two hits or overlapped by two hits.
        if args.keep_non_hits:
          if i == 0:
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = 1
            CutEnd = max(ThisEnd, NextStart)
          elif i == NumHits-1:
            LastStart, LastEnd = hits[i-1][5:7]
            CutStart = min(LastEnd, ThisStart) + 1
            CutEnd = SeqLength
          else:
            LastStart, LastEnd = hits[i-1][5:7]
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = min(LastEnd, ThisStart) + 1
            CutEnd = max(ThisEnd, NextStart)

        # For this block, we want to discard sequence that's not inside a hit,
        # and duplicate sequence that's overlapped by two hits.
        else:  
          CutStart = ThisStart
          CutEnd = ThisEnd

      else:

        # For this block, we want to keep sequence that's not inside a hit, and
        # cut in half sequence that's between two hits or overlapped by two
        # hits.
        if args.keep_non_hits:
          if i == 0:
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = 1
            CutEnd = (ThisEnd + NextStart) / 2
          elif i == NumHits-1:
            LastStart, LastEnd = hits[i-1][5:7]
            CutStart = (LastEnd + ThisStart) / 2 + 1
            CutEnd = SeqLength
          else:
            LastStart, LastEnd = hits[i-1][5:7]
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = (LastEnd + ThisStart) / 2 + 1
            CutEnd = (ThisEnd + NextStart) / 2

        # For this block, we want to discard sequence that's not inside a hit,
        # and cut in half sequence that's overlapped by two hits.
        else:
          if i == 0:
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = ThisStart
            CutEnd = min(ThisEnd, (ThisEnd + NextStart) / 2)
          elif i == NumHits-1:
            LastStart, LastEnd = hits[i-1][5:7]
            CutStart = max(ThisStart, (LastEnd + ThisStart) / 2 + 1)
            CutEnd = ThisEnd
          else:
            LastStart, LastEnd = hits[i-1][5:7]
            NextStart, NextEnd = hits[i+1][5:7]
            CutStart = max(ThisStart, (LastEnd + ThisStart) / 2 + 1)
            CutEnd = min(ThisEnd, (ThisEnd + NextStart) / 2)

      ThisCutSeq = copy.deepcopy(seq)
      ThisCutSeq.seq = ThisCutSeq.seq[CutStart-1 : CutEnd]
      ThisCutSeq.description = ''
      
      # Do we need to reverse complement?
      sstart, send = hit[7:9]
      if sstart > send:
        ThisCutSeq.seq = ThisCutSeq.seq.reverse_complement()

      ThisCutSeq.id += '.' + str(i + 1) + '_BlastsTo_' + \
      str(min(sstart, send)) + '-' + str(max(sstart, send))

      OutSeqs.append(ThisCutSeq)

# Write output
try:
  SeqIO.write(OutSeqs, args.out_file, "fasta")
except IOError:
  print('Problem writing output to', args.out_file)
  raise
