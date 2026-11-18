# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **282** | **280** | 2,3 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(282, 280), maxEE = c(2, 3), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 282 --p-trunc-len-r 280 \
    --p-max-ee-f 2 --p-max-ee-r 3
```

## What to expect

- About **70.3%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 550 nt. These settings can merge amplicons up to 550 nt, so your whole length range is covered.
- maxEE is 2,3, not DADA2's default of 2: at the recommended cut the 75th percentile of the reads' expected errors is 0.2 forward and 2.4 reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2). Only 64% of the read pairs pass at maxEE 2 and 70% pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable, because DADA2's error model uses the quality scores. The rule is AmpliTrim's own and has not been tested on independent data (see 'How maxEE is chosen' in the user guide).

## Things to check

- Read lengths vary (forward 70-289 nt, reverse 241-284 nt): DADA2 drops reads shorter than truncLen, which loses 10% of pairs here.
- Quality filtering removes 30% of read pairs at these settings.

## What maxEE keeps

At the recommended truncation (282 + 280), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 50.6% | 50.5% | 0.06 / 0.35 | - |
| 2 | 63.9% | 63.9% | 0.09 / 0.57 | - |
| **2,3 (used)** | 70.4% | 70.3% | 0.11 / 0.75 | 0.30 / 2.45 |
| 2,4 | 74.2% | 74.1% | 0.12 / 0.89 | 0.35 / 2.83 |
| 3 | 70.6% | 70.5% | 0.11 / 0.75 | 0.38 / 2.43 |
| 2,5 | 76.7% | 76.5% | 0.14 / 1.00 | 0.39 / 3.14 |
| 3,5 | 77.3% | 77.0% | 0.16 / 1.02 | 0.48 / 3.15 |
| 5 | 77.6% | 77.3% | 0.17 / 1.03 | 0.55 / 3.15 |

At maxEE 2,3, **70.4%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 49,114 pairs analysed (up to 289 + 284 nt).
- Amplicons: median 374 nt; 99% are 372-550 nt (all: 356-551), measured from the 91% of pairs that merge at full read length.
- Overlap at the recommended settings: about 187 nt on average.

![summary](summary.png)

## How these settings were chosen

- 49,114 of 589,361 read pairs were analysed. Totals below 562 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 562 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2,3) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.4.0_
