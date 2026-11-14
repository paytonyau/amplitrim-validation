# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **224** | **218** | 2 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(224, 218), maxEE = c(2, 2), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 224 --p-trunc-len-r 218 \
    --p-max-ee-f 2 --p-max-ee-r 2
```

## What to expect

- About **84.3%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 430 nt. These settings can merge amplicons up to 430 nt, so your whole length range is covered.
- maxEE is 2, DADA2's default: at the recommended cut the 75th percentile of the reads' expected errors is 0.4 forward and 0.4 reverse, both under 2, so the default already covers three quarters of the reads in each direction; 85% of the read pairs pass.

## Things to check

- Read lengths vary (forward 220-231 nt, reverse 209-228 nt): DADA2 drops reads shorter than truncLen, which loses 8% of pairs here.

## What maxEE keeps

At the recommended truncation (224 + 218), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 78.7% | 78.2% | 0.22 / 0.20 | - |
| **2 (used)** | 85.1% | 84.3% | 0.27 / 0.27 | - |
| 2,3 | 86.4% | 85.5% | 0.28 / 0.30 | 0.91 / 2.40 |
| 2,4 | 86.9% | 85.9% | 0.28 / 0.32 | 0.92 / 2.71 |
| 3 | 87.6% | 86.6% | 0.31 / 0.32 | 1.65 / 2.00 |
| 2,5 | 87.2% | 86.1% | 0.29 / 0.34 | 0.92 / 2.91 |
| 3,5 | 88.8% | 87.5% | 0.33 / 0.37 | 1.57 / 2.55 |
| 5 | 90.0% | 88.4% | 0.37 / 0.40 | 2.15 / 2.55 |

At maxEE 2, **85.1%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 43,981 pairs analysed (up to 231 + 228 nt).
- Amplicons: median 429 nt; 99% are 403-430 nt (all: 282-431), measured from the 97% of pairs that merge at full read length.
- Quality scores take only 4 distinct values, as with binned scores on newer Illumina instruments. DADA2's error model can need adjusting for such data.
- Overlap at the recommended settings: about 22 nt on average.

![summary](summary.png)

## How these settings were chosen

- 43,981 of 175,922 read pairs were analysed. Totals below 442 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 442 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
