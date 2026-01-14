# PCR Deduper

## Overview
PCR duplicates may remain in sequencing reads even after quality filtering, trimming, and alignment to a reference genome or transcriptome. The goal is to identify PCR duplicates and only keep one copy per read to reduce the source of noise in the data and ensure even amplification of DNA copies. 

This is a **Reference Based PCR Duplication Remover**. After the alignment step in a stranded sequencing workflow, the SAM file of uniquely mapped and single-ended reads can be sorted and checked for PCR duplicates using this tool.

 Given a SAM file sorted by chromosome and leftmost start position of uniquely mapped reads (that have been aligned already to a reference) that are single-ended, along with a set of known unique molecular indices (UMIs), this tool removes all PCR duplicates and only retains one copy per read, outputting a deduplicated SAM file.

## What is considered a duplicate?
For two reads to be considered duplicates that come from the same exact DNA or RNA molecule in our sample, the following conditions must be met:
- The reads are on the same chromosome
- The reads have the same 5' start position (after adjusting for soft-clipping and reference-consuming operations)
- The strands are the same
- The reads have matching UMIs
- Each read must have known UMIs (match a UMI in the set of known UMIs)

## Repository Structure
- `wang_deduper.py`: Python script that deduplicates single-ended reads.
- `unit_test/`: Directory containing sample data for developing and testing the pipeline, the expected output, and the documentation explaining all test cases.
- `test.sam`: Additional sorted SAM file used for developing and testing the pipeline.
- `STL96.txt`: Text file containing the list of 96 known valid UMIs.
- `pseudocode.md`: Draft and outline of the program logic.
- `report.txt`: Example report generated after running the script on `C1_SE_uniqAlign.sorted.sam`. 
Each run produces a report summarizing:
    - number of header lines
    - number of unique reads
    - number of incorrect UMIs
    - number of duplicates removed
    - number of reads per chromosome

## Requirements
- Python
- `argparse`, `re`

## Running Deduper
Use the help option to view all arguments and usage details:
`./wang_deduper.py -h`

The following command is an example of how to run this script:
`/usr/bin/time -v ./wang_deduper.py -f [sorted SAM file] -o [name of file to output] -u [text file containing all UMIs] -r [optional: text file containing reported values from the deduplicated SAM file]`