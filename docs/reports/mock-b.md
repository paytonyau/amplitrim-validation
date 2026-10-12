# Mock community B: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** 16S, V1-V3, MiSeq 2x300
- **Samples:** 8
- **Known truth:** Zymo reference sequences (three species, 7 variants; independent)
- **How we used it:** Independent test: kept aside and not used in the design
- **Data:** [PRJEB67827](https://www.ebi.ac.uk/ena/browser/view/PRJEB67827)
- **Source study:** reference [7] in the findings (to be verified)
- **Share of pairs passing DADA2's default filter (maxEE 2):** 20%

## In short

The poorest reads in the set: only 20% of the pairs pass DADA2's default filter. A looser maxEE gave a higher share of reads that match an expected sequence exactly (68.7% at maxEE 2, 73.6% to 75.2% for AmpliTrim's looser settings, 75.2% at maxEE 4,8 which is the released rule's setting, and 78.3% for the other tool's maxEE 6,9, which also found one more variant: 7 against 6). This is the one dataset where AmpliTrim's cut differed from the other tool's by more than the 2-point margin, by +2.6 points. In the first round of tests, which used a stricter margin of 1.0 point, it was also beyond the margin in AmpliTrim's favour.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 8 | +0.1 | -2.2 (-6.7 to +1.1) | 1.00 | +2.9 |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

*The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.*

| Samples | The other tool's cut and maxEE | AmpliTrim's cut | Difference in pairs merged (points) |
|---|---|---|---|
| 8 | 279/248, 6,9 | 274/253 | +2.6 |

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the rule changes the limit to **maxEE 4,8**. DADA2 merged this share of the read pairs:

*The share of the input read pairs that DADA2 merged with three maxEE settings, from the final build of the software.*

| Default maxEE 2 | maxEE 2,5 | Released rule | The rule's maxEE |
|---|---|---|---|
| 19.9% | 49.7% | 70.4% | 4,8 |

The rule has no cap. With a cap of 5 it would have set **maxEE 4,5** here, and DADA2 merged 55.0% of the pairs; without one it sets maxEE 4,8 and DADA2 merged 70.4%.

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +29.8 | noisier | more accurate |

Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:

| AmpliTrim setting | Difference from the other tool's own choice (points) |
|---|---|
| Default maxEE 2 | -54.1 |
| maxEE 2,5 | -24.2 |

**Accuracy on this mock community**

This dataset has expected sequences, so each setting can be scored.

| Setting | Variants found | Strains found | Reads matching an expected sequence exactly |
|---|---|---|---|
| AmpliTrim, maxEE 2 | 6 | 3 of 3 | 68.7% |
| AmpliTrim, maxEE 3,5 | 6 | 3 of 3 | 74.2% |
| AmpliTrim, maxEE 5 | 6 | 3 of 3 | 73.6% |
| AmpliTrim, maxEE 4,6 | 6 | 3 of 3 | 74.9% |
| The other tool, maxEE 6,9 | 7 | 3 of 3 | 78.3% |
| No truncation, maxEE 2 | 4 | 3 of 3 | 54.3% |

## Amplicon length and the longest amplicons

**The amplicon length the other tool needs**

The other tool needs you to tell it the amplicon length. We tried the longest (515 nt) and the typical (459 nt); it chose different cuts. AmpliTrim measures the lengths from your reads, so it needs no such number.

*The share of the input read pairs DADA2 merged with the cut each tool chose.*

| Length given: longest / typical | Other tool, given the longest | Other tool, given the typical | AmpliTrim, default |
|---|---|---|---|
| 515 / 459 | 279/248, 74.3% | 272/199, 45.7% | 19.9% |

**Leaving out the longest amplicons**

AmpliTrim sets the shortest cut that can merge your amplicons, but it leaves out the very longest few so that a handful of outliers do not force long, noisy reads. By default it leaves out the longest 0.2%. DADA2's merged share, at maxEE 2, when the longest 0.2%, 1% and 5% are left out:

*The share of the input read pairs DADA2 merged.*

| Longest 0.2% left out (default) | Longest 1% | Longest 5% |
|---|---|---|
| 19.9% | 19.9% | 19.9% |

## Which tests it took part in

This dataset was not part of: the test of choosing again on other samples.

## Where the numbers come from

Tables 1, 3, 4, 5, 6, 7, 8, 11 of the findings and the table of the final build on the validation page.

The web page of this dataset is in the results site of the archive: `datasets/mock-b/`.
