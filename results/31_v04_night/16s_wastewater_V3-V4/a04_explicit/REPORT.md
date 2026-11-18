# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **260** | **182** | 4,3 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(260, 182), maxEE = c(4, 3), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 260 --p-trunc-len-r 182 \
    --p-max-ee-f 4 --p-max-ee-r 3
```

## What to expect

- About **68.1%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 428 nt. These settings can merge amplicons up to 430 nt, so your whole length range is covered.
- maxEE is 4,3, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 3.3 forward and 2.1 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 55% of the read pairs pass at maxEE 2 and 75% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Only 2% of read pairs could be merged at full read length (amplicons longer than about 563 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 68% of read pairs are expected to pass filtering and merge.

## What maxEE keeps

At the recommended truncation (260 + 182), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 27.3% | 26.9% | 0.52 / 0.45 | - |
| 2 | 52.9% | 50.3% | 0.91 / 0.64 | - |
| 2,3 | 54.3% | 51.4% | 0.92 / 0.69 | 1.33 / 2.39 |
| 2,4 | 54.6% | 51.6% | 0.92 / 0.71 | 1.34 / 2.59 |
| 3 | 67.3% | 62.2% | 1.21 / 0.77 | 2.35 / 1.25 |
| 2,5 | 54.7% | 51.7% | 0.92 / 0.71 | 1.33 / 2.69 |
| **4,3 (used)** | 75.1% | 68.1% | 1.45 / 0.84 | 2.73 / 1.29 |
| 3,5 | 68.5% | 62.9% | 1.23 / 0.82 | 2.32 / 1.43 |
| 5 | 83.0% | 73.1% | 1.69 / 1.01 | 3.06 / 1.64 |

At maxEE 4,3, **75.1%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 49,440 pairs analysed (up to 291 + 284 nt).
- Amplicons: median 402 nt; 99% are 17-428 nt (all: 15-437), measured from the 2% of pairs that merge at full read length.
- Overlap at the recommended settings: about 28 nt on average.

![summary](summary.png)

## How these settings were chosen

- 49,440 of 2,224,789 read pairs were analysed. Totals below 440 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 440 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 4,3) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.4.0_
