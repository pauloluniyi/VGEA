#!/usr/bin/python
from __future__ import print_function

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script 'cleans' a sequence alignment. Currently it
uses input files required to be formatted in a very specific way, and so is
unlikely to be useful to anyone but the code's author, but I may make it more
flexible at some point. The steps are as follows. Any lower-case
base or "?" character is replaced by "N". The amplicon_regions_file given as an
argument groups the columns positions of the alignment (nominally into 'amplicon
regions'); the seq_based_blacklist given as an argument is used for specifying
which regions which for which sequences should be blacklisted (masked by a run
of "N"s). Additionally, any region that is more than 20% "N" characters is
blacklisted. The patient_based_blacklist given as an arg is used to wholly
remove the sequence(s) from particular individuals (nominally patients).
Multiple sequences for the same patient are flattened into one sequence by using
the base of the longest sequence (the most bases excluding N or gaps) if is not
an N, otherwise the base of the second-longest sequence if it is not an N, etc.
Then any gap character neighbouring an N is iteratively replaced by an N. Any
wholly undetermined (purely N) sequences are removed.'''

import argparse
import os
import sys
import itertools
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
import collections
from re import sub
import numpy as np
import matplotlib.pyplot as plt
from ShiverFuncs import TranslateSeqCoordsToAlnCoords

# Define a function to check files exist, as a type for the argparse.
def File(_file):
  if not os.path.isfile(_file):
    raise argparse.ArgumentTypeError(_file +' does not exist or is not a file.')
  return _file


# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('--global_aln', type=File, help='''An alignment of all
consensuses to be cleaned; use either this option or
It should contain no ambiguity codes; Ambiguity codes can be estimated using
shiver/tools/EstimateAmbiguousBases.py. If sequence names contain the string
"_consensus" and/or "_MinCov", the name is truncated by removing the string and
every thing after it. After that, each sequence name should be either the ID of
a patient, or the ID of a patient followed by an underscore then a unique
string, to distinguish multiple sequences associated with the same patient.''')
parser.add_argument('--individual_consensus', type=File, nargs='+',
help='''Use this option to specify one or more files containing an alignment of
a consensus to be cleaned and a reference sequence, with the same reference
sequence in all of these files. Use either this option or --global_aln. If you
use this option you must also use --common_reference.''')
parser.add_argument('--common_reference', help='''Use this to specify the name
of the reference sequence that is present in each of the individual consensus
files, if you are using the --individual_consensus option.''')
parser.add_argument('patient_based_blacklist', type=File, help='''A plain text
file in which each line is the name of a patient to be blacklisted, i.e. the
sequence(s) for that patient should be removed. A header line is not expected.
''')
parser.add_argument('seq_based_blacklist', type=File, help='''A csv-format file
specifying particular sequences and particular regions in sequences to be
blacklisted. A header is expected with the following column names: BAM (values
of this column should coincide with sequence names in the alignment, after the
latter have been truncated as described above), keep.at.all (values of this
column should be TRUE or FALSE, reflecting whether that sequence should be
wholly blacklisted), then should come the names of the different amplicon
regions (values of each of these columns should TRUE or FALSE, reflecting
whether that region for that sequence should be blacklisted), finally the column
'origin' is ignored.''')
parser.add_argument('amplicon_regions_file', type=File, help='''A csv file in
which the first column is the name of a region, the second column is the
(1-based) start position of that region, and the third column is the (1-based)
end position of that region. A header line is not expected. Regions should be
mutually exclusive and collectively exhaustive, i.e. every position in the
alignment should be in exactly one region.''')
parser.add_argument('--output', help='''Use this to specify where the output
will go. If you are using --global_aln, this should be a single file; if you are
using --individual_consensus, this should be a directory (in which we will
create one file per input individual consensus).''')
parser.add_argument('-V', '--verbose', action='store_true',
help="Print a line for each blacklisting action.")
parser.add_argument('-SA', '--split-amplicons', action='store_true',
help='''With this option we write a separate alignment for each region. Output
file names will have things appended to them appropriately. Note that patients
who have some but not all regions blacklisted will appear in the output for some
regions but not others.''')
parser.add_argument('--print_new_seq_blacklist', help='''Use this to specify the
name of file in which we'll store additions to the seq_based_blacklist (due to
our analysis here of completeness/missingness in each region).''')
parser.add_argument('--sanger_to_beehive_id_dict', type=File, help='''A csv file
with a header line; the first column is the sanger ID, the second column is the
BEEHIVE ID. If this option is used, input sequences are assumed to be labelled
by Sanger ID and their BEEHIVE ID will be looked up.''')
args = parser.parse_args()

# Check that we either have a global alignment or individual consensuses
have_global_aln = args.global_aln != None
have_individual_consensus = args.individual_consensus != None
if have_global_aln and have_individual_consensus:
  print('Use either --global_aln or --individual_consensus, not both.',
  'Quitting.', file=sys.stderr)
  exit(1)
if (not have_global_aln) and (not have_individual_consensus):
  print('You must use one of --global_aln or --individual_consensus.',
  'Quitting.', file=sys.stderr)
  exit(1)

# If we have individual consensus files, (a) the output arg should be a 
# directory - make it if it doesn't exist already; (b) --common_reference should
# have been specified.
if have_individual_consensus:
  if not os.path.isdir(args.output):
    os.mkdir(args.output)
  if args.common_reference == None:
    print('The --individual_consensus option requires the --common_reference',
    'option to be used too. Quitting.', file=sys.stderr)
    exit(1)

# Read in the sanger to BEEHIVE ID dict if desired.
have_sanger_dict = args.sanger_to_beehive_id_dict != None
if have_sanger_dict:
  sanger_to_beehive_id_dict = {}
  all_beehive_ids = set([])
  with open(args.sanger_to_beehive_id_dict, 'r') as f:
    for line_num_min_1, line in enumerate(f):
      if line_num_min_1 == 0:
        continue
      line = line.strip()
      fields = line.split(",")
      if len(fields) != 2:
        print('Encountered', len(fields), 'fields in', \
        args.sanger_to_beehive_id_dict + '. Expected 3. Quitting.',
        file=sys.stderr)
        exit(1)
      sanger_id, beehive_id = fields
      if sanger_id in sanger_to_beehive_id_dict:
        print('Encountered', sanger_id, 'a second time in', \
        args.sanger_to_beehive_id_dict + '. Quitting.', file=sys.stderr)
        exit(1)
      if beehive_id in all_beehive_ids:
        print('Encountered', beehive_id, 'a second time in', \
        args.sanger_to_beehive_id_dict + '. Quitting.', file=sys.stderr)
        exit(1)
      sanger_to_beehive_id_dict[sanger_id] = beehive_id
      all_beehive_ids.add(beehive_id)


# Read and check the amplicon regions. 
regions_dict = collections.OrderedDict()
regions_dict_not_empty = False
with open(args.amplicon_regions_file, 'r') as f:
  for line in f:
    fields = line.split(',')
    if len(fields) != 3:
      print('Encountered', len(fields), 'fields in', \
      args.amplicon_regions_file + '. Expected 3. Quitting.', file=sys.stderr)
      exit(1)
    region, start, end = fields
    start, end = int(start), int(end)
    if start > end:
      print('Start greater than end for region', region, 'in', 
      args.amplicon_regions_file + '. Quitting.', file=sys.stderr)
      exit(1)
    if region in regions_dict:
      print('Encountered', region, 'a second time in',
      args.amplicon_regions_file + '. Quitting.', file=sys.stderr)
      exit(1)
    if regions_dict_not_empty:
      previous_region = next(reversed(regions_dict))
      previous_start, previous_end = regions_dict[previous_region]
      if start != previous_end + 1:
        print('Region', region, 'starts at', start, 'which is not 1 more than',
        'the end of the previous region (' + str(previous_end), 'for',
        previous_region + '). Quitting.', file=sys.stderr)
        exit(1)
    #elif start != 1:
    #  print('The first region in', args.amplicon_regions_file,
    #  'does not start at 1. Quitting.', file=sys.stderr)
    #  exit(1)
    regions_dict[region] = (start, end)
    regions_dict_not_empty = True
last_region = next(reversed(regions_dict))
last_start, last_end = regions_dict[last_region]
regions = regions_dict.keys()
num_regions = len(regions)

# If we're splitting amplicons, set up file names for the per-region output
# files.
if args.split_amplicons:

  # If an output fasta file was specified, use its extension one was specified
  # (splitting into two pieces at the right-most dot).
  # If an output directory was specified, we'll choose the full file names
  # within that directory.
  # Check file names don't exist already.
  if have_global_aln:
    pieces = args.output.rsplit(".", 1)
    if len(pieces) == 2:
      per_region_output_file_dict = {region : pieces[0] + "_" + region + "." + \
      pieces[1] for region in regions}
    else:
      per_region_output_file_dict = {region : args.output + "_" + region + \
      ".fasta" for region in regions}
  else:
    per_region_output_file_dict = {region : os.path.join(args.output,
    "region_" + region + ".fasta") for region in regions}
  for output_file in per_region_output_file_dict.values():
    if os.path.isfile(output_file):
      print(output_file, "exists already; quitting to prevent overwriting.",
      file=sys.stderr)
      exit(1)
      

# Read the seq-based blacklist.
seq_blacklist_dict = {}
expected_header_line = 'BAM,keep.at.all,' + ','.join(regions) + ',origin'
num_fields = expected_header_line.count(',') + 1
with open(args.seq_based_blacklist, 'r') as f:
  for line_num_min_1, line in enumerate(f):

    # Check for the expected header
    if line_num_min_1 == 0:
      line = line.strip()
      if line != expected_header_line:
        print('Unexpected header line\n' + line + '\nfor',
        args.seq_based_blacklist + '; expected\n' + \
        expected_header_line + '\nQuitting.', file=sys.stderr)
        exit(1)
      continue

    # Check for the right number of fields
    fields = line.split(',')
    if len(fields) != num_fields:
      print('Unexpected number of fields', len(fields), 'on line',
      line_num_min_1 + 1, 'in', args.seq_based_blacklist + '; expected',
      str(num_fields) + '. Quitting.', file=sys.stderr)
      exit(1)

    # Coerce string bools to bools, and record them.
    values = fields[1:-1]
    if any(value != "TRUE" and value != "FALSE" for value in values):
      print('Unexpected value not equal to TRUE or FALSE on line',
      line_num_min_1 + 1, 'in', args.seq_based_blacklist + '. Quitting.',
      file=sys.stderr)
      exit(1)
    values = np.array([True if value == "TRUE" else False for value in values],
    dtype=bool)

    # If we've seen this seq in the blacklist already, we should blacklist each
    # region if either entry says so, i.e. only keep it if they agree that we
    # should.
    id_ = fields[0]
    if id_ in seq_blacklist_dict:
      previous_values = seq_blacklist_dict[id_]
      for i, new_value in enumerate(values):
        previous_value = previous_values[i]
        if new_value and previous_value:
          values[i] = True
        else:
          values[i] = False

    seq_blacklist_dict[id_] = values

blacklisted_patients = set([])
with open(args.patient_based_blacklist, 'r') as f:
  for line in f:
    patient = line.strip()
    blacklisted_patients.add(patient)

def preprocess_seq(seq, check_for_ambig_bases):
  '''Replaces lower-case bases and '?' by N, optionally checks for amibuous
  bases, checks whether it's wholly undetermined.'''
  seq_as_str = str(seq.seq)
  initial_length = len(seq_as_str)
  seq_as_str = sub("[a-z]|\?", "N", seq_as_str)
  assert len(seq_as_str) == initial_length
  if check_for_ambig_bases and any(not base in "ACGTN-" for base in seq_as_str):
    print('Seq', seq.id, 'contains a base other than A, C, G, T, N or -;',
    'unexpected. Have you run shiver/tools/EstimateAmbiguousBases.py on your',
    'alignment first? Quitting.', file=sys.stderr)
    exit(1)
  wholly_undetermined = all(base == "N" or base == "-" for base in seq_as_str)
  return seq_as_str, wholly_undetermined

def preprocess_seq_id(seq_id):
  "Removes things we don't want from the end of seq ids."
  seq_id = seq_id.split('_consensus', 1)[0]
  seq_id = seq_id.split('_MinCov', 1)[0]
  seq_id = seq_id.split('_remap', 1)[0]
  return seq_id

def get_beehive_id(seq_id):
  "Return whatever's before the first underscore if there is one."
  return seq_id.split("_", 1)[0]

def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Iteratively replace gaps that border N by N.

  Where N neighbours a gap, propagate N outwards until it touches bases on both
  sides (because deletions should only be called when the bases on either side
  are known). e.g.
  ACTG---N---ACTG
  becomes
  ACTGNNNNNNNACTG'''
  
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == 'N':
      BaseToLeftIsNoCoverage = True
      ResultingSeq += 'N'
    elif base == '-':
      if BaseToLeftIsNoCoverage:
        ResultingSeq += 'N'
      else:
        ResultingSeq += '-'
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  assert not "N-" in ResultingSeq and not "-N" in ResultingSeq
  return ResultingSeq

def num_known_bases(seq):
  return len(seq) - seq.count("N") - seq.count("-")

if have_global_aln:

  # Read the alignment
  try:
    alignment = AlignIO.read(args.global_aln, "fasta")
  except:
    print('Problem reading', args.global_aln + ':', file=sys.stderr)
    raise
  alignment_length = alignment.get_alignment_length()

  # Check that the end of the last region is the length of the alignment.
  #if last_end != alignment_length:
  #  print('The last region in', args.amplicon_regions_file,
  #  'does not end at the alignment length (' + str(alignment_length) + \
  #  '). Quitting.', file=sys.stderr)
  #  exit(1)

  # collection_of_seqs is the thing we'll iterate through for seq processing.
  # For the global alignment, that's the seq objects therein; for individual
  # consensus files, that the set of files.
  collection_of_seqs = alignment

  # With a global alignment we'll check for ambiguous bases (because these can
  # be much more easily estimated before hand using
  # ~/shiver/tools/EstimateAmbiguousBases.py) and record a dict of seqs as we
  # iterate through them.
  check_ambig_bases = True
  seq_dict = collections.OrderedDict()

else:
  collection_of_seqs = args.individual_consensus
  check_ambig_bases = False
  seq_dict = set([]) # Used just to record what seqs we've seen so far.

extra_blacklisted_seq_ids = set([])
first_gapless_reference = None
blacklisted_patients_found = set([])
blacklisted_seqs_found = set([])
max_missingness = 0.2
unaln_seqs_by_region = {region : {} for region in regions}
#completeness_percents = collections.defaultdict(list)

for seq in collection_of_seqs:

  # If individual consensus files were given, seq is one such file. Read the
  # file and reassign seq to be the first sequence in it. 
  if have_individual_consensus:
    individual_consensus_file = seq
    try:
      individual_alignment = AlignIO.read(individual_consensus_file, "fasta")
    except:
      print('Problem reading', individual_consensus_file + ':', file=sys.stderr)
      raise
    alignment_length = individual_alignment.get_alignment_length()
    seq = individual_alignment[0]

    # Get the reference out of the individual consensus file. Translate the
    # region coords with respect to the reference to be with respect to this
    # alignment. Check the reference has the same sequence (after removing gaps)
    # in every file, overall and for each separate region.
    this_reference = None
    for seq_ in individual_alignment:
      if seq_.id == args.common_reference:
        this_reference = str(seq_.seq)
        this_reference_gapless = this_reference.replace("-", "")
    if this_reference == None:
      print('Could not find', args.common_reference, 'in',
      individual_consensus_file + '. Quitting.', file=sys.stderr)
      exit(1)
    regions_dict_this_aln = {}
    this_reference_by_region = {}
    for region, region_boundaries in regions_dict.items():
      start_in_this_aln, end_in_this_aln = \
      TranslateSeqCoordsToAlnCoords(this_reference, region_boundaries)
      regions_dict_this_aln[region] = (start_in_this_aln, end_in_this_aln)
      this_reference_by_region[region] = \
      this_reference[start_in_this_aln - 1 : end_in_this_aln].replace("-", "") 
    if first_gapless_reference == None:
      first_gapless_reference = this_reference_gapless
      first_reference_by_region = this_reference_by_region
    else:
      if this_reference_gapless != first_gapless_reference:
        print('The sequence of', args.common_reference, 'differs between',
        args.individual_consensus[0], 'and', individual_consensus_file + \
        '. Quitting.', file=sys.stderr)
        exit(1)
      for region in regions:
        if this_reference_by_region[region] != \
        first_reference_by_region[region]:
          print('Error: internal malfunction of', __file__ + ': problem',
          'translating coordinates. the sequence of', args.common_reference,
          'in region', region, 'differs between', args.individual_consensus[0],
          'and', individual_consensus_file + '. Quitting.', file=sys.stderr)
          exit(1)

  else:
    # If we're working with the global alignment, the coordinates of the regions
    # are the same for every seq.
    regions_dict_this_aln = regions_dict

  # Preprocess the seq id and check we haven't seen it before.
  seq_id = preprocess_seq_id(seq.id)
  if seq_id in seq_dict:
    print('Encountered seq', seq_id, 'a second time. Quitting.',
    file=sys.stderr)
    exit(1)

  # Look up the beehive_sanger style ID from the sanger ID.
  if have_sanger_dict:
    if seq_id in sanger_to_beehive_id_dict:
      seq_id = sanger_to_beehive_id_dict[seq_id]
    else:
      print('Error: seq', seq_id, 'not found in the first column of', + \
      args.sanger_to_beehive_id_dict + '. Quitting.', file=sys.stderr)
      exit(1)

  # Preprocess the seq, and skip if it's wholly undetermined.
  seq_as_str, wholly_undetermined = preprocess_seq(seq, check_ambig_bases)
  if wholly_undetermined:
    if args.verbose:
      print("Skipping sequence", seq_id, "which is wholly undetermined after",
      "removing lower-case bases.")
    continue

  # Add our own blacklisting of regions, based on the fraction that's not "N",
  # for regions that have not already been blacklisted.
  for region_num in xrange(num_regions):
    seq_in_blacklist = seq_id in seq_blacklist_dict
    if seq_in_blacklist and not seq_blacklist_dict[seq_id][region_num + 1]:
      continue
    region = regions[region_num]
    start, end = regions_dict_this_aln[region]
    region_length = end - start + 1
    seq_here = seq_as_str[start - 1: end]
    missingness = float(seq_here.count("N")) / region_length
    if missingness >= max_missingness:
      if not seq_in_blacklist:
        seq_blacklist_dict[seq_id] = [True] * (num_regions + 1)
      seq_blacklist_dict[seq_id][region_num + 1] = False
      extra_blacklisted_seq_ids.add(seq_id)
    #completeness_percents[region].append(completeness_percent)

  # Delete every seq from a blacklisted patient      
  beehive_id = get_beehive_id(seq_id)
  if beehive_id in blacklisted_patients:
    if args.verbose:
      print("Sequence ", seq_id, ": discarding whole sequence as patient ",
      beehive_id, " was blacklisted.", sep='')
    blacklisted_patients_found.add(beehive_id)
    continue

  if seq_id in seq_blacklist_dict:

    blacklisted_seqs_found.add(seq_id)

    # The seq_blacklist_values are bools saying keep this seq at all, then keep
    # each region.
    seq_blacklist_values = seq_blacklist_dict[seq_id]
    keep_at_all = seq_blacklist_values[0]
    if not keep_at_all:
      if args.verbose:
        print("Sequence ", seq_id,
        ": discarding whole sequence as it was blacklisted.", sep='')
      continue

    for region_num in xrange(num_regions):

      # Mask blacklisted regions by "N".
      # Remember start and end are 1-based coords.
      keep_region = seq_blacklist_values[region_num + 1]
      if not keep_region:
        region = regions[region_num]
        start, end = regions_dict_this_aln[region]
        if args.verbose:
          print("Sequence ", seq_id, ": masking region ", region + \
          " (positions ", start, "-", end, ").", sep='')
        seq_as_str = seq_as_str[:start - 1] + \
        "N" * (end - start + 1) + seq_as_str[end:]

  # If we have individual consensuses, either record each region of this
  # consensus into a separate list for later aggregation between patients if 
  # desired, otherwise just write the output file.
  if have_individual_consensus:
    seq_as_str = PropagateNoCoverageChar(seq_as_str)
    if args.split_amplicons:
      for region, (start, end) in regions_dict_this_aln.items():
        region_length = end - start + 1
        seq_here = seq_as_str[start - 1: end].replace("-", "")
        if not all(base == "N" for base in seq_here):
          if beehive_id in unaln_seqs_by_region[region]:
            previous_seq_here = unaln_seqs_by_region[region][beehive_id]
            if num_known_bases(seq_here) <= num_known_bases(previous_seq_here):
              if args.verbose:
                print("Ignoring", seq_id, "in region", region, "because a",
                "previously encountered seq from", beehive_id, "has at least",
                "as many known bases.")
              continue
            if args.verbose:
              print("Using", seq_id, "in region", region, "instead of a",
              "previously encountered seq from", beehive_id, "because the",
              "former has more known bases.")
          unaln_seqs_by_region[region][beehive_id] = seq_here
    if all(base == "N" for base in seq_as_str):
      if args.verbose:
        print("Skipping sequence", seq_id, "which is wholly undetermined after",
        "blacklisting.")
      continue
    out_file = os.path.join(args.output, seq_id + ".fasta")
    if os.path.isfile(out_file):
      print(out_file, "exists already; quitting to prevent overwriting.",
      file=sys.stderr)
      exit(1)
    individual_alignment[0].seq = Seq.Seq(seq_as_str)
    individual_alignment[0].id = seq_id
    individual_alignment[0].description = ""
    AlignIO.write(individual_alignment, out_file, "fasta")
    seq_dict.add(seq_id)

  # If we have a global alignment, just record the seq for later processing.
  else:
    seq_dict[seq_id] = seq_as_str

# Warn about blacklisted patients that were not found.
blacklisted_patients_not_found = \
blacklisted_patients - blacklisted_patients_found
if len(blacklisted_patients_not_found) != 0:
  print('Warning: for the following blacklisted patients, specified in',
  args.patient_based_blacklist + ', no sequence was found:',
  ' '.join(blacklisted_patients_not_found), file=sys.stderr)

# Warn about blacklisted seqs that were not found.
blacklisted_seqs_not_found = \
set(seq_blacklist_dict.keys()) - blacklisted_seqs_found
if len(blacklisted_seqs_not_found) != 0:
  print('Warning: the following blacklisted seqs, specified in',
  args.seq_based_blacklist + ', were not encountered:',
  ' '.join(blacklisted_seqs_not_found), file=sys.stderr)
  
# If desired, write the new blacklist as a result of our missingness analysis.
if args.print_new_seq_blacklist != None:
  with open(args.print_new_seq_blacklist, 'w') as f:
    f.write(expected_header_line.rsplit(',', 1)[0] + "\n")
    for seq_id in sorted(extra_blacklisted_seq_ids):
      seq_blacklist_values = seq_blacklist_dict[seq_id]
      f.write(seq_id + "," + ",".join(map(str, seq_blacklist_values)).upper() \
      + "\n")

#for region, completeness_percents_here in completeness_percents.items():
#  fig, ax = plt.subplots()
#  ax.set_yscale("log", nonposy='clip')
#  ax.set_ylim(bottom=0.5, top=len(alignment))
#  plt.xlabel("Percent of alignment positions that are not 'N'", fontsize=15)
#  plt.ylabel("Number of sequences", fontsize=15)
#  plt.hist(completeness_percents_here, bins=100)
#  plt.savefig("AmpliconCompletenessHisto_" + region + ".pdf")
#  plt.clf()
#exit(0)

# There's nothing more to do for individual consensus files, except generating
# region-specific files if desired.
if have_individual_consensus:
  if args.split_amplicons:
    for region, output_file in per_region_output_file_dict.items():
      seq_dict_here = unaln_seqs_by_region[region]
      sorted_seqs_here = sorted(seq_dict_here.items(), key=lambda x:x[0])
      SeqIO.write((SeqIO.SeqRecord(Seq.Seq(seq), id=seq_id, description='') \
      for seq_id, seq in sorted_seqs_here), output_file, 'fasta')
  exit(0)

# Form a dict for patients with multiple sequences, listing the sequence IDs
per_patient_seq_counts = \
collections.Counter(get_beehive_id(seq_id) for seq_id in seq_dict)
multiply_seqd_patients = set(patient for patient, count in \
per_patient_seq_counts.items() if count > 1)
singly_seqd_patients = set(patient for patient, count in \
per_patient_seq_counts.items() if count == 1)
seq_ids_for_multiply_seqd_patients = collections.defaultdict(list)
for seq_id in seq_dict:
  beehive_id = get_beehive_id(seq_id)
  if beehive_id in multiply_seqd_patients:
    seq_ids_for_multiply_seqd_patients[beehive_id].append(seq_id)


# Merge multiple seqs per patient.
for multiply_seqd_patient, seq_ids in \
seq_ids_for_multiply_seqd_patients.items():
  seqs = [seq_dict[seq_id] for seq_id in seq_ids]
  num_seqs = len(seqs)
  seqs_by_length = sorted(seqs, key=lambda seq : alignment_length - \
  seq.count("N") - seq.count("-"), reverse=True)
  best_seq = ""
  for pos in xrange(alignment_length):

    # Try each seq in order from best to worst until one of them has something
    # other than an N:
    best_base = "N"
    for i in xrange(num_seqs):
      base = seqs_by_length[i][pos]
      if base != "N":
        best_base = base
        break
    best_seq += best_base

  for seq_id in seq_ids:
    del seq_dict[seq_id]
  seq_dict[multiply_seqd_patient] = best_seq

# Multiply sequenced patients now have only one seq, labelled by the patient;
# ensure that for other patients their only seq is also labelled by the patient.
seq_ids_to_rename = \
[seq_id for seq_id in seq_dict if get_beehive_id(seq_id) != seq_id]
for seq_id in seq_ids_to_rename:
  seq_dict[get_beehive_id(seq_id)] = seq_dict.pop(seq_id)

# Iteratively replace gaps that border N by N.
for seq_id in seq_dict:
  seq_dict[seq_id] = PropagateNoCoverageChar(seq_dict[seq_id])

# Sort seqs by name
sorted_seqs = sorted(seq_dict.items(), key=lambda x:x[0])

if args.split_amplicons:

  for region, (start, end) in regions_dict.items():
    OutSeqs = []
    region_length = end - start + 1
    for seq_id, seq in sorted_seqs:
      seq_here = seq[start - 1: end]
      if not all(base == "N" for base in seq_here):
        OutSeqs.append(SeqIO.SeqRecord(Seq.Seq(seq_here), id=seq_id,
        description=''))
    SeqIO.write(OutSeqs, per_region_output_file_dict[region], 'fasta')    

OutSeqs = []
for seq_id, seq in sorted_seqs:
  if all(base == "N" for base in seq):
    if args.verbose:
      print("Skipping sequence", seq_id, "which is wholly undetermined after",
      "blacklisting.")
  else:
    OutSeqs.append(SeqIO.SeqRecord(Seq.Seq(seq), id=seq_id, description=''))
SeqIO.write(OutSeqs, args.output, 'fasta')
        

