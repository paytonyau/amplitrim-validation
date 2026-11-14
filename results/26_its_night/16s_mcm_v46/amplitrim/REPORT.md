# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **273** | **265** | 4,3 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(273, 265), maxEE = c(4, 3), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 273 --p-trunc-len-r 265 \
    --p-max-ee-f 4 --p-max-ee-r 3
```

## What to expect

- About **60.8%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 526 nt. These settings can merge amplicons up to 526 nt, so your whole length range is covered.
- maxEE is 4,3, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 3.0 forward and 2.6 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 43% of the read pairs pass at maxEE 2 and 65% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Only 89% of read pairs could be merged at full read length (amplicons longer than about 548 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 61% of read pairs are expected to pass filtering and merge.

## What maxEE keeps

At the recommended truncation (273 + 265), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 19.0% | 18.8% | 0.37 / 0.37 | - |
| 2 | 43.2% | 42.1% | 0.63 / 0.63 | - |
| 2,3 | 51.5% | 49.4% | 0.58 / 0.93 | 0.34 / 2.46 |
| 2,4 | 56.6% | 53.2% | 0.57 / 1.15 | 0.40 / 2.84 |
| 3 | 59.6% | 56.7% | 0.84 / 0.85 | 1.39 / 1.42 |
| 2,5 | 59.6% | 55.1% | 0.58 / 1.32 | 0.45 / 3.14 |
| **4,3 (used)** | 65.0% | 60.8% | 1.05 / 0.81 | 1.89 / 1.18 |
| 3,5 | 68.2% | 62.6% | 0.82 / 1.23 | 1.14 / 2.25 |
| 5 | 78.0% | 69.4% | 1.20 / 1.15 | 1.92 / 1.80 |

At maxEE 4,3, **65.0%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 24,899 pairs analysed (up to 279 + 281 nt).
- Amplicons: median 523 nt; 99% are 521-526 nt (all: 283-527), measured from the 89% of pairs that merge at full read length.
- Overlap at the recommended settings: about 15 nt on average.

![summary](summary.png)

## How these settings were chosen

- 24,899 of 24,899 read pairs were analysed. Totals below 538 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 538 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 4,3) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
