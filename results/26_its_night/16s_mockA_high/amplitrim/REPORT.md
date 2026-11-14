# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **273** | **247** | 2,3 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(273, 247), maxEE = c(2, 3), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 273 --p-trunc-len-r 247 \
    --p-max-ee-f 2 --p-max-ee-r 3
```

## What to expect

- About **71.6%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 508 nt. These settings can merge amplicons up to 508 nt, so your whole length range is covered.
- maxEE is 2,3, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 1.6 forward and 2.4 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 65% of the read pairs pass at maxEE 2 and 72% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Only 89% of read pairs could be merged at full read length (amplicons longer than about 582 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Quality filtering removes 28% of read pairs at these settings.

## What maxEE keeps

At the recommended truncation (273 + 247), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 43.4% | 43.4% | 0.20 / 0.43 | - |
| 2 | 63.8% | 63.8% | 0.35 / 0.72 | - |
| **2,3 (used)** | 71.7% | 71.6% | 0.40 / 0.91 | 0.81 / 2.44 |
| 2,4 | 75.5% | 75.4% | 0.43 / 1.03 | 0.86 / 2.77 |
| 3 | 75.4% | 75.2% | 0.50 / 0.94 | 1.33 / 2.19 |
| 2,5 | 77.6% | 77.4% | 0.45 / 1.13 | 0.88 / 3.03 |
| 3,5 | 83.1% | 82.8% | 0.58 / 1.21 | 1.32 / 2.85 |
| 5 | 87.8% | 87.2% | 0.76 / 1.29 | 1.82 / 2.82 |

At maxEE 2,3, **71.7%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 47,380 pairs analysed (up to 297 + 297 nt).
- Amplicons: median 484 nt; 99% are 460-505 nt (all: 378-536), measured from the 89% of pairs that merge at full read length.
- Overlap at the recommended settings: about 41 nt on average.

![summary](summary.png)

## How these settings were chosen

- 47,380 of 142,140 read pairs were analysed. Totals below 520 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 520 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2,3) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
