# Definitions used in the tables and the register

Every table of recall in this archive (sections 33 to 45, the paper, the workbook and the web pages) uses these words in this way.

| Word | Meaning |
|---|---|
| **Reference sequence** | A known sequence of the mock community: an isolate's or a strain's Sanger sequence, or an NCBI reference for the ATCC standard. |
| **Member** | A reference sequence, or a group of reference sequences within 2 edits of each other (they cannot be told apart). 189 Pauvert strains give 151 expected members; Bakker's two *Fusarium* isolates that differ by one base are one member. Partial references (Pauvert) are judged as prefixes. |
| **True ASV** | An ASV within 2 edits of a reference sequence. |
| **False ASV, false reads** | An ASV farther than 2 edits from every reference sequence, and the reads in it. |
| **Found** | A member is found in a sample when a true ASV matches it and holds at least 0.1% of the sample's reads (the rule of sections 33 to 44). Section 45 also counts a member as found with any read at all. |
| **Recall** | The share of the members expected in a sample that are found in it. |
| **Yield, reads kept** | The share of the input read pairs left after DADA2's filtering, denoising, merging and chimera removal. |

**What "found" is not.** It is a match of sequences, not a taxonomic assignment: no reference database (UNITE or other) is used and no rank (kingdom to species) is involved. An ASV that is correct to genus or species but more than 2 edits from every reference counts as false, so the measure is stricter than a genus or species call. Strains with identical sequences cannot be told apart, so a count of members is not a count of species. Pauvert's 76.5% of strains (against the authors' 75.8%) credits every strain that shares a found member's sequence; the counting rules of the two studies are not guaranteed identical.

For the 16S mock communities (section 45) a strain is found when one of its variants is found exactly, as `mock_check.py` judges.
