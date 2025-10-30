#!/usr/bin/env python

# make sure the input SAM file is sorted...

# while writing the code, 
# run with: ./wang_deduper.py -u STL96.txt -f unit_test/sorted_input.sam -o temp_unit_test_output.sam

# remember to bin usr/time -v when actually running!

#--------------
# import modules
#--------------

import argparse
import re

#--------------
# argparse: ADD REQUIRED = TRUE LATER
#--------------

def get_args():
    parser = argparse.ArgumentParser(description="to deduplicate single-ended reads given a SAM file of uniquely mapped reads and a file with a list of UMIs")
    parser.add_argument("-f", help="designates file path to sorted SAM file", type=str)
    parser.add_argument("-o", help="designates file path to deduplicated SAM file", type=str)
    parser.add_argument("-u", help="designates file containing the list of UMIs", type=str)
    parser.add_argument("-r", help="designated file path to report values from the deduplicated SAME file", type=str)
    return parser.parse_args()

args = get_args()

# Set global argparse variables:
file = args.f
outfile = args.o
umi = args.u

#--------------
# functions
#--------------

def strandedness(SAM_col_2):

    '''
    Given the bitwise FLAG (SAM col 2) in the SAM file, return FALSE if the strand is the forward 
    or return TRUE if the strand is reverse.
    '''

    SAM_col_2 = int(SAM_col_2) # making sure it's an integer value for bitwise comparison

    # TRUE indicates a reverse strand
    if ((SAM_col_2 & 16) == 16):
        rev_comp = True
    else:
        # FALSE indicates a forward strand
        rev_comp = False

    return(rev_comp)

def extract_read_info(SAM_line):
    '''
    Given a SAM line (already stripped and split by tabs), extract the chromosome, 
    strand (using strandedness(), which returns FALSE for a forward strand or TRUE for a reverse strand), 
    non-adjusted start position, and cigar string. 
    Return them as separate values.
    '''

    # get chr: col 3
    chrom = SAM_line[2]

    # get strand: col 2
    strand = strandedness(SAM_line[1])  # FALSE for a forward strand, TRUE for a reverse strand

    # get nonadj_pos, aka the 1-based leftmost position: col 4
    nonadj_pos = int(SAM_line[3])

    # get cigar, a string: col 6
    cigar = SAM_line[5]

    return chrom, strand, nonadj_pos, cigar

def get_adj_pos(this_strand, nonadj_pos, cigar) -> int:

    '''
    Takes:
    whether the strand is a forward or reverse strand (this_strand),
    the "nonadjusted" 1-based, leftmost starting position of the read (SAM col 4) (nonadj_pos),
    and the cigar string (cigar).

    Checks the cigar string to see if soft-clipping (S) is occuring,
    and to see if there are instances of the reference being consumed
    (such as deletions (D), skipped regions (N), mismatches/matches (M)).

    Returns:
    The calculated actual, adjusted 5' position based on the following rules:
    - if is a Forward strand:
        - subtract leading soft clipping (S at the start of the cigar string)
    - if is a Reverse strand:
        - add deletions (D) and skipped regions (N) and mismatches/matches (M)
        - DO NOT add insertions (I), which do not consume reference
        - then add trailing soft clipping (S at the end of the cigar string)
    '''
    
    # parse through cigar string

    cigar_tuples = re.findall()
    # no valid CIGAR info
    if not cigar_tuples =
        return nonadj_pos

    if this_strand == FALSE: # forward strand
        # if there is soft clipping, extract the soft clip-- make sure it's NOT the #S at the end of the cigar string!
        soft_clip_value = <>
        # adj_pos = nonadj_pos - soft_clip_value
    
    if this_strand == TRUE: # reverse strand
        # if there is soft clipping, extract the soft clip-- make sure it IS the #S at the end of the cigar string!
        soft_clip_value = <>

        # consider things in the cigar that consumes reference!!!
        #M, the total matched/unmatched bases
        strand_len = <>
        #D, the deletions
        deletion = <>
        #N, the gaps
        gaps = <>

        # adj_pos = nonadj_pos + strand_len + + deletion + gaps + soft_clip_value

    return adj_pos

#--------------
# initialize things
#--------------

# Initialize a set for known UMIs
UMIs = set()

# Read in the known UMIs text file and add each line to the set of known UMIs
# umi argparse variable is set to STL96.txt
with open(umi, 'r') as f:
    for line in f:
        line = line.strip()
        UMIs.add(line) # add each line to the set of known UMIs

#--------------
# code
#--------------

# Open a file for writing the output SAM file
with open(outfile, "w") as o:
    # Read in the sorted SAM file
    with open(file, 'r') as f:
        # Track which chromosome I'm on to parse through the reads by chromosome
        current_chr = None
        # Initialize a set for storing reads I've parsed through already: (UMI, chrom, strand, pos)
        seen_read = set() # This should reset per chromosome

        # Parse through each line of input SAM file
        for line in f:
            # if line is a header line:
            if line.startswith("@"):
                # write line to output file
                o.write(line)
                continue # skip the rest if the line is a header line

            # if line is not a header line:
            line = line.strip().split("\t") # strip and split based on tab separation
            this_umi = line[0].split(":")[7]
            if this_umi not in UMIs:
                continue # skip the rest if UMI is not in the set of known UMIs

            # get this_chr, this_strand, nonadj_pos, and cigar-- from two functions that I wrote
            this_chr, this_strand, nonadj_pos, cigar = extract_read_info(line)

            if this_chr != current_chr:
                seen_read = set() # reset seen reads if this_chr is "new"
                current_chr = this_chr # update current_chr to the "new" chromosome

            # get the actual, adjusted, 5' starting pos using a function I wrote
            this_pos = get_adj_pos(this_strand, nonadj_pos, cigar)

            # add all these to temporary storage for this read
            read_id = (this_umi, this_chr, this_strand, this_pos)

            # identify if this is a PCR duplicate!
            if read_id not in seen_read:
                # if read_id is not in seen_read, write it to the output
                seen_read.add(read_id)
                # write the current line as a read in the output file, tab separated and with new lines
                o.write("\t".join(line) + "\n")
            else:
                # skip PCR duplicate!
                continue

# regex... +