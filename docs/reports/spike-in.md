# Mouse gut with spike-ins, V3-V4: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** 16S, V3-V4, NovaSeq 6000, 2x250, binned qualities; mouse gut with synthetic spike-ins
- **Samples:** 12
- **Known truth:** spike-in sequences not available, so not used
- **How we used it:** Independent test: kept aside and not used in the design
- **Data:** [PRJNA1243306](https://www.ebi.ac.uk/ena/browser/view/PRJNA1243306)
- **Source study:** to add
- **Share of pairs passing DADA2's default filter (maxEE 2):** 88%

## In short

Mouse gut with synthetic spike-ins; the spike-in sequences were not available, so accuracy was not assessed. NovaSeq reads with binned quality scores. At the same maxEE, AmpliTrim's cut merged 1.3 points more than the other tool's (+0.4 to +2.4). In the first round of tests, which used a stricter margin of 1.0 point, this was beyond the margin in AmpliTrim's favour.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 12 | -0.5 | +3.4 (+2.9 to +3.7) | 1.00 | +3.4 |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

*The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.*

| Samples | The other tool's cut and maxEE | AmpliTrim's cut | Difference in pairs merged (points) |
|---|---|---|---|
| 12 | 232/207, 2,1 | 230/210 | +1.3 (+0.4 to +2.4) |

*The numbers in brackets are a 95% interval from resampling the samples; it is given only when there are 10 or more samples.*

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the released rule keeps DADA2's default of 2, as the validation page records for the clean datasets and for ITS.

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +2.3 | noisier | no expected sequences |

Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:

| AmpliTrim setting | Difference from the other tool's own choice (points) |
|---|---|
| Default maxEE 2 | +6.3 |
| maxEE 2,5 | +8.6 |

## Would the choice change on other samples?

Each tool chose a cut again on three separate sets of samples from this dataset. The choice counts as the same if the cuts differ by at most 10 nt and the maxEE is the same. Cuts are forward/reverse; the maxEE is in brackets.

| AmpliTrim's cuts | Same each time? | The other tool's cuts | Same each time? |
|---|---|---|---|
| 222/218; 231/208; 229/210 (2) | yes | 222/217 (1); 210/229 (1,2); 233/206 (2,1) | no |

## Which tests it took part in

This dataset was not part of: the accuracy scoring (it needs expected sequences); the test of the amplicon length given to the other tool; the test of leaving out the longest amplicons.

## Where the numbers come from

Tables 1, 3, 4, 5, 9, 11 of the findings.

The web page of this dataset is in the results site of the archive: `datasets/spike-in/`.
