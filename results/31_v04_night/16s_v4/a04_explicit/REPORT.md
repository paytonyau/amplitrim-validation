# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **124** | **115** | 2 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(124, 115), maxEE = c(2, 2), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 124 --p-trunc-len-r 115 \
    --p-max-ee-f 2 --p-max-ee-r 2
```

## What to expect

- About **98.5%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 227 nt. These settings can merge amplicons up to 227 nt, so your whole length range is covered.
- maxEE is 2, DADA2's default: at the recommended cut the 75th percentile of the reads' expected errors is 0.0 forward and 0.1 reverse, both under 2, so the default already covers three quarters of the reads in each direction; 99% of the read pairs pass.

## What maxEE keeps

At the recommended truncation (124 + 115), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 96.5% | 96.3% | 0.07 / 0.08 | - |
| **2 (used)** | 98.7% | 98.5% | 0.09 / 0.10 | - |
| 2,3 | 98.9% | 98.7% | 0.09 / 0.10 | 0.75 / 2.34 |
| 2,4 | 99.0% | 98.8% | 0.09 / 0.10 | 0.64 / 2.64 |
| 3 | 99.3% | 99.0% | 0.10 / 0.10 | 1.74 / 1.74 |
| 2,5 | 99.1% | 98.8% | 0.09 / 0.11 | 0.54 / 3.06 |
| 3,5 | 99.5% | 99.1% | 0.10 / 0.11 | 1.41 / 2.31 |
| 5 | 99.8% | 99.2% | 0.11 / 0.12 | 2.16 / 2.42 |

At maxEE 2, **98.7%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 49,720 pairs analysed (up to 230 + 248 nt).
- Amplicons: median 207 nt; 99% are 206-227 nt (all: 139-241), measured from the 97% of pairs that merge at full read length.
- Overlap at the recommended settings: about 32 nt on average.

![summary](summary.png)

## How these settings were chosen

- 49,720 of 3,032,869 read pairs were analysed. Totals below 239 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 239 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.4.0_
