# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **195** | **168** | 2 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(195, 168), maxEE = c(2, 2), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 195 --p-trunc-len-r 168 \
    --p-max-ee-f 2 --p-max-ee-r 2
```

## What to expect

- About **90.6%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 351 nt. These settings can merge amplicons up to 351 nt, so your whole length range is covered.
- maxEE is 2, DADA2's default: at the recommended cut the 75th percentile of the reads' expected errors is 0.2 forward and 0.3 reverse, both under 2, so the default already covers three quarters of the reads in each direction; 91% of the read pairs pass.

## Things to check

- Only 80% of read pairs could be merged at full read length (amplicons longer than about 582 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.

## What maxEE keeps

At the recommended truncation (195 + 168), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 84.5% | 84.4% | 0.11 / 0.14 | - |
| **2 (used)** | 90.7% | 90.6% | 0.16 / 0.20 | - |
| 2,3 | 92.0% | 91.8% | 0.18 / 0.24 | 0.94 / 2.42 |
| 2,4 | 92.6% | 92.3% | 0.18 / 0.25 | 0.93 / 2.71 |
| 3 | 93.9% | 93.6% | 0.22 / 0.26 | 1.84 / 1.85 |
| 2,5 | 92.8% | 92.5% | 0.18 / 0.27 | 0.89 / 2.95 |
| 3,5 | 95.0% | 94.4% | 0.23 / 0.30 | 1.69 / 2.36 |
| 5 | 97.2% | 96.2% | 0.32 / 0.35 | 2.45 / 2.34 |

At maxEE 2, **90.7%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 49,235 pairs analysed (up to 297 + 297 nt).
- Amplicons: median 315 nt; 99% are 297-346 nt (all: 98-547), measured from the 80% of pairs that merge at full read length.
- Overlap at the recommended settings: about 45 nt on average.

![summary](summary.png)

## How these settings were chosen

- 49,235 of 1,427,804 read pairs were analysed. Totals below 363 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 363 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.4.0_
