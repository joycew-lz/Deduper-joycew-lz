# Input file input.sam:
- Is sorted (by chrom (col 3) and leftmost mapping position (col 4)).

- Lines 1-24: headers
- Line 25, Line 26: PCR duplicates
- Line 27, Line 28: Not PCR duplicates: Different chromosomes
- Line 29, Line 30: Not PCR duplicates: Different strands
- Line 31, Line 32: PCR duplicates due to soft clipping-- note the different 1-based, leftmost starting positions, but the same adjusted 5' starting positions, for FORWARD strands (adj_pos = nonadj_pos - soft_clip_value, where soft_clip_value is NOT the #S at the end of the cigar string).
- Line 33, Line 34: PCR duplicates due to soft clipping-- note the different 1-based, leftmost starting positions, but the same adjusted 5' starting positions, for REVERSE strands (adj = nonadj_pos + strand_len + soft_clip_value, where soft_clip_value IS the #S at the end of the cigar string, and strand_len is the total number of matched or unmatched bases #M).
- Line 35, Line 36: Not PCR duplicates: Different UMIs
- Line 37, Line 38: Not PCR duplicates: Though they could be PCR duplicates due to soft clipping as FORWARD strands, Line 37 has a UMI that doesn't belong in the set of known UMIs.
- New line character at the end.

# Output file output.sam:
- Is sorted (by chrom (col 3) and leftmost mapping position (col 4)).

- Lines 1-24: headers written out
- Line 25: keep only Line 25 from input file: keep only 1 copy of the read when encountering PCR duplicates
- Line 26, Line 27: keep both Line 27 and Line 28 from input file
- Line 28, Line 29: keep both Line 29 and Line 30 from input file
- Line 30: keep only Line 31 from input file: keep only 1 copy of the read when encountering PCR duplicates due to soft clipping for FORWARD strands.
- Line 31: keep only Line 33 from input file: keep only 1 copy of the read when encountering PCR duplicates due to soft clipping for REVERSE strands.
- Line 32, Line 33: keep both Line 35 and Line 36 from input file
- Line 34: only keep Line 38 from input file, since Line 37 can be tossed, as it has a UMI that doesn't belong in the set of known UMIs.
- New line character at the end.