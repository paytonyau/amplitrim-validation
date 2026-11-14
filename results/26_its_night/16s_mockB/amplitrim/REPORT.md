# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **274** | **253** | 4,8 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(274, 253), maxEE = c(4, 8), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 274 --p-trunc-len-r 253 \
    --p-max-ee-f 4 --p-max-ee-r 8
```

## What to expect

- About **67.2%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 515 nt. These settings can merge amplicons up to 515 nt, so your whole length range is covered.
- maxEE is 4,8, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 3.5 forward and 7.1 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 20% of the read pairs pass at maxEE 2 and 72% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Only 50% of read pairs could be merged at full read length (amplicons longer than about 565 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 67% of read pairs are expected to pass filtering and merge.
- The rule sets maxEE 4,8, above 5: these reads are very noisy; more low-quality reads may not help. Check on a mock community if you can.

## What maxEE keeps

At the recommended truncation (274 + 253), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 5.3% | 5.3% | 0.22 / 0.66 | - |
| 2 | 18.3% | 18.2% | 0.35 / 1.26 | - |
| 2,3 | 31.9% | 31.7% | 0.42 / 1.79 | 0.50 / 2.49 |
| 2,4 | 43.1% | 42.7% | 0.47 / 2.23 | 0.56 / 2.94 |
| 3 | 33.0% | 32.8% | 0.48 / 1.80 | 0.64 / 2.46 |
| 2,5 | 50.9% | 50.1% | 0.52 / 2.57 | 0.61 / 3.30 |
| 3,5 | 54.1% | 53.2% | 0.63 / 2.62 | 0.78 / 3.31 |
| 5 | 56.1% | 55.0% | 0.74 / 2.64 | 0.93 / 3.31 |
| **4,8 (used)** | 72.2% | 67.2% | 0.91 / 3.48 | 1.10 / 4.23 |

At maxEE 4,8, **72.2%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 34,140 pairs analysed (up to 289 + 288 nt).
- Amplicons: median 465 nt; 99% are 431-515 nt (all: 319-522), measured from the 50% of pairs that merge at full read length.
- Overlap at the recommended settings: about 59 nt on average.

![summary](summary.png)

## How these settings were chosen

- 34,140 of 34,140 read pairs were analysed. Totals below 527 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 527 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 4,8) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
