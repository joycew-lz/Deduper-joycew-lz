#!/usr/bin/env python

# make sure the input SAM file is sorted.

# while writing the code, 
# run with: ./wang_deduper.py -u STL96.txt -f unit_test/sorted_input.sam -o temp_unit_test_output.sam

# remember to bin usr/time -v

#--------------
# import modules
#--------------

import argparse

#--------------
# argparse: ADD REQUIRED = TRUE LATER
#--------------

def get_args():
    parser = argparse.ArgumentParser(description="to deduplicate single-ended reads given a SAM file of uniquely mapped reads and a file with a list of UMIs")
    parser.add_argument("-f", help="designates absolute file path to sorted SAM file", type=str)
    parser.add_argument("-o", help="designates absolute file path to deduplicated SAM file", type=str)
    parser.add_argument("-u", help="designates file containing the list of UMIs", type=str)
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
    Given the bitwise FLAG (SAM col 2) in the SAM file, return FALSE if the strand is the forward or return TRUE if the strand is reverse.
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
    Given a SAM line (already stripped and split by tabs), extract the chromosome, strand (using strandedness(), which returns FALSE for a forward strand or TRUE for a reverse strand), non-adjusted start position, and cigar string. Return them as separate values.
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

# Initialize a set for storing reads I've parsed through already: (UMI, chrom, strand, pos)
# This should reset to be empty after I get through each chromosome in the SAM file!
seen_read = set()

# Open a file for writing the output SAM file
with open(outfile, "w") as o:
    # Read in the sorted SAM file
    with open(file, 'r') as f:
        # Parse through each line of input SAM file
        for line in f:
            # if line is a header line:
            if line.startswith("@"):
                # write line to output file
                o.write(line)

            # if line is not a header line:
            else:
                line = line.strip().split("\t") # strip and split based on tab separation
                this_umi = line[0].split(":")[7]
                if this_umi not in UMIs:
                    continue # skip if UMI is not in the set of known UMIs!

                # get this_chrom, this_strand, nonadj_pos, cigar from two functions that I wrote!
                this_chr, this_strand, nonadj_pos, cigar = extract_read_info(line)


# regex... +