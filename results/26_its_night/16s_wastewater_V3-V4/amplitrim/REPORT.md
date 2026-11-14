# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **261** | **183** | 4,3 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(261, 183), maxEE = c(4, 3), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 261 --p-trunc-len-r 183 \
    --p-max-ee-f 4 --p-max-ee-r 3
```

## What to expect

- About **65.6%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 428 nt. These settings can merge amplicons up to 432 nt, so your whole length range is covered.
- maxEE is 4,3, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 3.4 forward and 2.3 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 54% of the read pairs pass at maxEE 2 and 73% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Only 2% of read pairs could be merged at full read length (amplicons longer than about 561 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 66% of read pairs are expected to pass filtering and merge.

## What maxEE keeps

At the recommended truncation (261 + 183), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 23.8% | 23.6% | 0.54 / 0.45 | - |
| 2 | 49.9% | 47.6% | 0.95 / 0.68 | - |
| 2,3 | 51.5% | 48.8% | 0.96 / 0.73 | 1.33 / 2.37 |
| 2,4 | 51.9% | 49.1% | 0.97 / 0.75 | 1.34 / 2.57 |
| 3 | 64.8% | 59.9% | 1.27 / 0.82 | 2.34 / 1.30 |
| 2,5 | 52.0% | 49.2% | 0.97 / 0.76 | 1.34 / 2.70 |
| **4,3 (used)** | 72.5% | 65.6% | 1.50 / 0.88 | 2.72 / 1.34 |
| 3,5 | 66.2% | 60.8% | 1.29 / 0.88 | 2.31 / 1.50 |
| 5 | 80.7% | 70.5% | 1.75 / 1.06 | 3.04 / 1.69 |

At maxEE 4,3, **72.5%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 44,423 pairs analysed (up to 288 + 285 nt).
- Amplicons: median 402 nt; 99% are 17-428 nt (all: 13-428), measured from the 2% of pairs that merge at full read length.
- Overlap at the recommended settings: about 28 nt on average.

![summary](summary.png)

## How these settings were chosen

- 44,423 of 44,423 read pairs were analysed. Totals below 440 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 440 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 4,3) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
