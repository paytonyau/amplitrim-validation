# Mock community A: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** 16S, V1-V3 (27F-YM/534R), MiSeq 2x300
- **Samples:** 6 high-input samples for the settings (23 in the dilution series)
- **Known truth:** the authors' ASVs anchored to Sanger sequences of 14 of 15 strains (not independent)
- **How we used it:** Used in the design: used while designing AmpliTrim, so not independent evidence
- **Data:** [PRJEB67797](https://www.ebi.ac.uk/ena/browser/view/PRJEB67797)
- **Source study:** reference [7] in the findings (to be verified)
- **Share of pairs passing DADA2's default filter (maxEE 2):** 65%

## In short

We used this dataset while designing AmpliTrim, so it is not independent evidence. It has six high-input samples. At every maxEE setting we tried, the same variants and strains were recovered, so a looser maxEE cost no accuracy here. Skipping truncation altogether merged 23.3 points fewer pairs.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 6 | -0.8 | +2.7 (+0.3 to +4.2) | 0.90 | -0.4 |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

*The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.*

| Samples | The other tool's cut and maxEE | AmpliTrim's cut | Difference in pairs merged (points) |
|---|---|---|---|
| 6 | 273/247, 3,4 | 273/247 | 0.0 |

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the rule changes the limit to **maxEE 2,3**. DADA2 merged this share of the read pairs:

*The share of the input read pairs that DADA2 merged with three maxEE settings, from the final build of the software.*

| Default maxEE 2 | maxEE 2,5 | Released rule | The rule's maxEE |
|---|---|---|---|
| 62.7% | 76.1% | 69.5% | 2,3 |

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +13.4 | noisier | as accurate |

Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:

| AmpliTrim setting | Difference from the other tool's own choice (points) |
|---|---|
| Default maxEE 2 | -15.3 |
| maxEE 2,5 | -1.8 |

**Accuracy on this mock community**

This dataset has expected sequences, so each setting can be scored.

| Setting | Variants found | Strains found | Reads matching an expected sequence exactly |
|---|---|---|---|
| AmpliTrim at each maxEE tried, and the other tool | 36 | 12 of 14 | 99.7% to 99.8% |
| No truncation, maxEE 2 | 33 | 11 of 14 | 99.8% |

## Which tests it took part in

This dataset was not part of: the test of the amplicon length given to the other tool; the test of leaving out the longest amplicons; the test of choosing again on other samples.

## Where the numbers come from

Tables 1, 3, 4, 5, 6, 11 of the findings and the table of the final build on the validation page.

The web page of this dataset is in the results site of the archive: `datasets/mock-a/`.
