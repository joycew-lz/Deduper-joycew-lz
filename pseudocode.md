# Reference-based PCR Duplicate Removal Assignment

# Goal:
Given a SAM file of uniquely mapped reads (that have been aligned already to a reference) that are single-ended, along with a set of known UMIs, create a tool that will remove all PCR duplicates and only retain one copy per read.

# Define the problem:
PCR duplicates may exist in the reads even after the quality filtering/trimming and alignment to a reference genome/transcriptome steps of DNA or RNA sequencing. We want to identify PCR duplicates and only keep one copy per read to reduce the source of noise in our data and ensure even amplification of DNA copies. After the alignment step in a stranded sequencing workflow, we can take the SAM file of uniquely mapped and single-ended reads, sort the SAM file, and check for PCR duplicates.

For two reads to be considered duplicates that come from the same exact DNA or RNA molecule in our sample, we must check:
- That the reads are on the same chromosome
- That the reads have the same 5' start position (to adjust for soft clipping)
- That the strands are the same
- That the reads have matching UMIs, or unique molecular indices
    - (If the reads appear to not have a known UMI, don't consider retaining the read)

# Write examples:
Include a properly formated sorted input sam file
Include a properly formated expected output sam file

# Develop your algorithm using pseudocode

```
import argparse

# Argparse things

def get_args():
    parser = argparse.ArgumentParser(description="to deduplicate single-ended reads given a SAM file of uniquely mapped reads and a file with a list of UMIs")
    parser.add_argument("-f", help="designates absolute file path to sorted sam file", type=str, required=True)
    parser.add_argument("-o", help="designates absolute file path to deduplicated sam file", type=str, required=True)
    parser.add_argument("-u", help="designates file containing the list of UMIs", type=str, required=True)
    return parser.parse_args()

args = get_args()

# Set global argparse variables:
file = args.f
outfile = args.o
umi = args.u

# Initialize a set for known UMIs
UMIs = set()

# Read in the known UMIs text file and add each line to the set of known UMIs
# umi argparse variable is set to STL96.txt
with open(umi, 'r') as f:
    for line in f:
        line = line.strip()
        UMIs.add(line) # add each line to the set of known UMIs
    close umi

# Sort the SAM file using samtools sort

# Initialize a set for storing reads I've parsed through already: (UMI, chrom, strand, pos)
seen_read = set()

# Open a file for writing the output SAM file
with open(outfile, "w") as o:
    # Read in the sorted SAM file
    with open(file, 'r') as f:
        #  Parse through each line of input SAM file...
        for line in f:
            line = line.strip().split("\t") # strip and split based on tab separation
            # if the line is a header
                # write line to output file
            
            else: # if the line is not a header
                # check that the UMI is in the set of known UMIs
                this_umi = line[0]
                if this_umi not in UMIs:
                    continue # skip if UMI is not in the set of known UMIs!

                # get this_chrom, this_strand, nonadj_pos, cigar from two functions that I wrote!
                this_chr, this_strand, nonadj_pos, cigar = extract_read_info(line)

                # check to see if soft clipping is happening, and to get the actual, adjusted, 5' starting pos:
                this_pos = get_adj_pos(this_strand, nonadj_pos, cigar)
            
                # add all these to temporary storage
                current_read = (this_umi, this_chr, this_strand, this_pos)

                # identify if this is a PCR duplicate!
                if current_read not in seen_read:
                    seen_read.add(current_read)
                    o.write("\t".join(line) + "\n") # write the current line as the first read in the output file, tab separated and with new lines
                else:
                    continue # skip PCR duplicate!
        
```


# Determine high level functions

```
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

Input: 16
Output: True

Input: 83
Output: True

Input: 0
Output: False
```

```
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

Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU

Output: 2, False, 76814284, 71M
```

```
def get_adj_pos(this_strand, nonadj_pos, cigar) -> int:

    '''
    Take  whether the strand is a forward or reverse strand (this_strand) and the "nonadjusted" 1-based, leftmost starting position of the read (SAM col 4) (nonadj_pos)
    to see if there is soft clipping happening (SAM col 6, indicated by 'S'). 
    Calculate the actual, adjusted, 5' starting position.
    '''
    
    # check to see if soft clipping is present in the cigar string
    <>
        if not:
        adj_pos = nonadj_pos

    if this_strand == FALSE: # forward strand
        # if there is soft clipping, extract the soft clip-- make sure it's NOT the #S at the end of the cigar string!
            soft_clip_value = <>
        # adj_pos = nonadj_pos - soft_clip_value
    
    if this_strand == TRUE: # reverse strand
        # if there is soft clipping, extract the soft clip-- make sure it IS the #S at the end of the cigar string!
            soft_clip_value = <>
        strand_len = <> # #M, the total matched bases
        # adj_pos = nonadj_pos + strand_len + soft_clip_value

    return adj_pos

Input: FALSE, 10, 3S5M2S
Output: 7

Input: TRUE, 10, 3S5M2S
Output: 15
```