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
        #  Parse through each line of input SAM file...
        for line in f:
            if line.startswith("@"):
                o.write(line)
                # write line to output file
            else:
                line = line.strip().split("\t") # strip and split based on tab separation
                this_umi = line[0].split(":")[7]
                print(this_umi)