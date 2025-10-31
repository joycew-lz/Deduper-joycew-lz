#!/usr/bin/env python

# make sure the input SAM file is sorted...
# command line: cp /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam .
    # samtools sort C1_SE_uniqAlign.sam -o C1_SE_uniqAlign.sorted.sam
    # removed original file

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
    parser.add_argument("-f", help="designates file path to sorted SAM file", type=str, required=True)
    parser.add_argument("-o", help="designates file path to deduplicated SAM file", type=str, required=True)
    parser.add_argument("-u", help="designates file containing the list of UMIs", type=str, required=True)
    parser.add_argument("-r", help="designated file path to report values from the deduplicated SAM file", type=str, required=False)
    return parser.parse_args()

args = get_args()

# Set global argparse variables:
file = args.f
outfile = args.o
umi = args.u
report = args.r

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
    
    # Parse through cigar string first:
    cigar_list = re.findall(r'\d+[MIDNS]', cigar) # cigar_list output is a list of strings, like ['3S', '5M', '10N', '2S']

    # now, need to separate the number from the letter;
    # set up a list of tuples that will store integer(s)-string pairs, like [(3, 'S'), (5, 'M'), (10, 'N'), (2, 'S')]
    cigar_list_tuples = []
    # Parse through cigar_list;
        # for a number(s)-letter combo (that are both characters in a string), like '3S'
    for num_letter_combo in cigar_list:
        numbers = int(num_letter_combo[:-1])# extract everything but the last character (the number(s)), and make sure it's an integer
        letter = num_letter_combo[-1] # extract the last character (the letter)
        # and append it to the list of tuples, as a tuple (number(s), letter)
        cigar_list_tuples.append((numbers, letter))

    if this_strand == False: # forward strand
        if cigar_list_tuples[0][1] == "S":
            leading_soft = cigar_list_tuples[0][0]
        else:
            leading_soft = 0
        adj_pos = nonadj_pos - leading_soft
    
    if this_strand == True: # reverse strand
        reference_consuming = {"M", "D", "N"}
        reference_consuming_total = 0
        for number, letter in cigar_list_tuples:
            if letter in reference_consuming:
                reference_consuming_total += number
        if cigar_list_tuples[-1][1] == "S":
            trailing_soft = cigar_list_tuples[-1][0]
        else:
            trailing_soft = 0
        adj_pos = nonadj_pos + trailing_soft + reference_consuming_total

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

# Initialize counters
header_lines = 0
unique_reads = 0
wrong_UMIs = 0
duplicates_removed = 0
counts_per_chr = {}

#--------------
# code
#--------------

# Open a file for writing the output SAM file
with open(outfile, "w") as o:
    # Read in the sorted SAM file
    with open(file, 'r') as f:
        # Track which chromosome I'm on to parse through the reads by chromosome
        current_chr = ""
        # Initialize a set for storing reads I've parsed through already: (UMI, chrom, strand, pos)
        seen_reads = set() # This should reset per chromosome

        # Parse through each line of input SAM file
        for line in f:
            # if line is a header line:
            if line.startswith("@"):
                header_lines += 1
                # write line to output file
                o.write(line)
                continue # skip the rest if the line is a header line

            # if line is not a header line:
            line = line.strip().split("\t") # strip and split based on tab separation
            this_umi = line[0].split(":")[-1]
            if this_umi not in UMIs:
                wrong_UMIs += 1
                continue # skip the rest if UMI is not in the set of known UMIs

            # get this_chr, this_strand, nonadj_pos, and cigar-- from two functions that I wrote
            this_chr, this_strand, nonadj_pos, cigar = extract_read_info(line)

            if this_chr != current_chr:
                seen_reads = set() # reset seen reads if this_chr is "new"
                current_chr = this_chr # update current_chr to the "new" chromosome

            # get the actual, adjusted, 5' starting pos using a function I wrote
            this_pos = get_adj_pos(this_strand, nonadj_pos, cigar)

            # add all these to temporary storage for this read
            read_id = (this_umi, this_chr, this_strand, this_pos)

            # identify if this is a PCR duplicate!

            # if read_id isn't a duplicate:
            if read_id not in seen_reads:
                unique_reads += 1
                # if read_id is not in seen_reads, write it to the output
                seen_reads.add(read_id)

                # add to counts per chromosome
                if this_chr not in counts_per_chr:
                    counts_per_chr[this_chr] = 0
                counts_per_chr[this_chr] += 1

                # write the current line as a read in the output file, tab separated and with new lines
                o.write("\t".join(line) + "\n")
            
            # skip PCR duplicate!
            else:
                duplicates_removed += 1
                continue

#--------------
# report
#--------------

if args.r:
    with open(report, "w") as r:
        r.write(
            f"Deduplicating for: {file}\n"
            f"Number of header lines: {header_lines}\n"
            f"Number of unique reads: {unique_reads}\n"
            f"Number of wrong UMIs: {wrong_UMIs}\n"
            f"Number of duplicates removed: {duplicates_removed}\n"
        )
        r.write(f"\nReads per chromosome:\n")
        for chrom, count in sorted(counts_per_chr.items()):
            r.write(f"{chrom}\t{count}\n")

# test with: ./wang_deduper.py -f unit_test/sorted_input.sam -o unit_test/temp_unit_test_actual_output.sam -u STL96.txt
# run with: 
    # /usr/bin/time -v ./wang_deduper.py -f C1_SE_uniqAlign.sorted.sam -o C1_SE_uniqAlign.sorted_output.sam -u STL96.txt -r report.txt
