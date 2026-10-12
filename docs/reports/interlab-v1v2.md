# Inter-laboratory mock, V1-V2: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** 16S, V1-V2 (about 337 bp), MiSeq 2x300
- **Samples:** 6
- **Known truth:** none (15 bacteria; digital PCR family abundances only)
- **How we used it:** Independent test: kept aside and not used in the design
- **Data:** [PRJEB34919](https://www.ebi.ac.uk/ena/browser/view/PRJEB34919)
- **Source study:** reference [10] in the findings (to be verified)
- **Share of pairs passing DADA2's default filter (maxEE 2):** 96%

## In short

Six samples with no known truth (15 bacteria; only digital PCR abundances of families). The reads are clean: 96% of the pairs pass DADA2's default filter. When the other tool was given the typical amplicon length instead of the longest, it merged 54.4% of the pairs against 93.6%.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 6 | -0.1 | +0.3 (-0.1 to +0.4) | 1.00 | +0.2 |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

*The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.*

| Samples | The other tool's cut and maxEE | AmpliTrim's cut | Difference in pairs merged (points) |
|---|---|---|---|
| 6 | 156/169, 1 | 155/170 | 0.0 |

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the released rule keeps DADA2's default of 2, as the validation page records for the clean datasets and for ITS.

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +1.4 | noisier | no expected sequences |

Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:

| AmpliTrim setting | Difference from the other tool's own choice (points) |
|---|---|
| Default maxEE 2 | +2.3 |
| maxEE 2,5 | +3.7 |

## Amplicon length and the longest amplicons

**The amplicon length the other tool needs**

The other tool needs you to tell it the amplicon length. We tried the longest (313 nt) and the typical (295 nt); it chose different cuts. AmpliTrim measures the lengths from your reads, so it needs no such number.

*The share of the input read pairs DADA2 merged with the cut each tool chose.*

| Length given: longest / typical | Other tool, given the longest | Other tool, given the typical | AmpliTrim, default |
|---|---|---|---|
| 313 / 295 | 156/169, 93.6% | 147/160, 54.4% | 95.9% |

**Leaving out the longest amplicons**

AmpliTrim sets the shortest cut that can merge your amplicons, but it leaves out the very longest few so that a handful of outliers do not force long, noisy reads. By default it leaves out the longest 0.2%. DADA2's merged share, at maxEE 2, when the longest 0.2%, 1% and 5% are left out:

*The share of the input read pairs DADA2 merged.*

| Longest 0.2% left out (default) | Longest 1% | Longest 5% |
|---|---|---|
| 95.9% | 95.9% | 93.8% |

## Which tests it took part in

This dataset was not part of: the accuracy scoring (it needs expected sequences); the test of choosing again on other samples.

## Where the numbers come from

Tables 1, 3, 4, 5, 7, 8, 11 of the findings.

The web page of this dataset is in the results site of the archive: `datasets/interlab-v1v2/`.
