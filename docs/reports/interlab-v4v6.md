# Inter-laboratory mock, V4-V6: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** 16S, V4-V6 (about 564 bp), MiSeq 2x300
- **Samples:** 5 of 6
- **Known truth:** none
- **How we used it:** Independent test: kept aside and not used in the design (real tight overlap)
- **Data:** [PRJEB34919](https://www.ebi.ac.uk/ena/browser/view/PRJEB34919)
- **Source study:** reference [10] in the findings (to be verified)
- **Share of pairs passing DADA2's default filter (maxEE 2):** 43%

## In short

A real case of tight overlap: amplicons of about 564 bp. Only 43% of the pairs pass DADA2's default filter, and relaxing to maxEE 2,5 merged 18.8 points more.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 5 | -0.2 | -1.6 (-3.8 to +0.4) | 1.00 | +0.9 |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

*The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.*

| Samples | The other tool's cut and maxEE | AmpliTrim's cut | Difference in pairs merged (points) |
|---|---|---|---|
| 5 | 270/268, 4 | 271/267 | +0.3 |

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the rule changes the limit to **maxEE 4,3**. DADA2 merged this share of the read pairs:

*The share of the input read pairs that DADA2 merged with three maxEE settings, from the final build of the software.*

| Default maxEE 2 | maxEE 2,5 | Released rule | The rule's maxEE |
|---|---|---|---|
| 41.6% | 60.4% | 63.4% | 4,3 |

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +18.8 | mixed | no expected sequences |

Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:

| AmpliTrim setting | Difference from the other tool's own choice (points) |
|---|---|
| Default maxEE 2 | -27.0 |
| maxEE 2,5 | -8.2 |

## Which tests it took part in

This dataset was not part of: the accuracy scoring (it needs expected sequences); the test of the amplicon length given to the other tool; the test of leaving out the longest amplicons; the test of choosing again on other samples.

## Where the numbers come from

Tables 1, 3, 4, 5, 11 of the findings and the table of the final build on the validation page.

The web page of this dataset is in the results site of the archive: `datasets/interlab-v4v6/`.
