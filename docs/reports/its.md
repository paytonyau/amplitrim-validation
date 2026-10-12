# Spring wheat, ITS: evidence report

*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*

- **Marker, region, platform:** ITS (ITS1-Fl2/ITS2), MiSeq 2x250
- **Samples:** 37 spring-wheat samples
- **Known truth:** none
- **How we used it:** ITS case study: a case study of the --its option (the other tool cannot be run)
- **Data:** [PRJEB58189](https://www.ebi.ac.uk/ena/browser/view/PRJEB58189)
- **Source study:** reference [8] in the findings (to be verified)
- **Share of pairs passing DADA2's default filter (maxEE 2):** 74%

## In short

The ITS case study (the --its option). The other tool cannot be run, because it needs reads of one length. The filter estimate was within 0.1 points. The merge estimate was +21.1 points too high (with 5,000 pairs per sample). Raising maxEE to 3 gained 3.7 points on 37 of 37 samples.

## How good was AmpliTrim's estimate?

Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.

*How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.*

| Samples | Pass-rate estimate off by (points) | Merge estimate off by (points; average and range over the settings tried) | Same order as DADA2 (1.00 = identical) | The other tool's pass-rate estimate off by (points) |
|---|---|---|---|---|
| 37 | +0.1 | +21.1 (+13.7 to +30.4) | 1.00 | cannot be run |

The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.

## Compared with the other tool, at the same maxEE

The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.

The other tool cannot be run on this dataset, because it needs reads of one length.

## Choosing the maxEE limit

maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.

**What the released rule does**

On this dataset the released rule keeps DADA2's default of 2, as the validation page records for the clean datasets and for ITS.

**What a looser limit would do**

We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.

*Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.*

| Change in pairs merged (points) | Quality of the extra reads kept | Accuracy (needs expected sequences) |
|---|---|---|
| +2.7 | not known | no expected sequences |

## Which tests it took part in

This dataset was not part of: the comparison with the other tool at the same maxEE; the accuracy scoring (it needs expected sequences); the test of the amplicon length given to the other tool; the test of leaving out the longest amplicons; the test of choosing again on other samples.

## Where the numbers come from

Tables 1, 3, 11 of the findings.

The web page of this dataset is in the results site of the archive: `datasets/its/`.
