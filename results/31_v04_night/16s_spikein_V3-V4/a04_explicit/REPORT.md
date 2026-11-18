# AmpliTrim report

## Recommended settings

| truncLenF | truncLenR | maxEE | minOverlap |
|---|---|---|---|
| **229** | **210** | 2 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR,
                     truncLen = c(229, 210), maxEE = c(2, 2), truncQ = 2, maxN = 0)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 229 --p-trunc-len-r 210 \
    --p-max-ee-f 2 --p-max-ee-r 2
```

## What to expect

- About **88.5%** of read pairs should pass quality filtering and merge.
- Your longest amplicon is 427 nt. These settings can merge amplicons up to 427 nt, so your whole length range is covered.
- maxEE is 2, DADA2's default: at the recommended cut the 75th percentile of the reads' expected errors is 0.7 forward and 0.5 reverse, both under 2, so the default already covers three quarters of the reads in each direction; 90% of the read pairs pass.

## What maxEE keeps

At the recommended truncation (229 + 210), the share of read pairs that pass DADA2's quality filter, and that then pass filtering and merge, for other maxEE values (forward, reverse):

| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|
| 1 | 79.5% | 78.8% | 0.33 / 0.25 | - |
| **2 (used)** | 89.6% | 88.5% | 0.43 / 0.33 | - |
| 2,3 | 91.0% | 89.7% | 0.44 / 0.36 | 1.21 / 2.40 |
| 2,4 | 91.4% | 90.0% | 0.45 / 0.37 | 1.22 / 2.64 |
| 3 | 93.9% | 92.3% | 0.51 / 0.39 | 2.04 / 1.67 |
| 2,5 | 91.6% | 90.1% | 0.45 / 0.38 | 1.21 / 2.78 |
| 3,5 | 95.2% | 93.2% | 0.53 / 0.44 | 2.01 / 2.15 |
| 5 | 97.9% | 95.1% | 0.61 / 0.49 | 2.60 / 2.23 |

At maxEE 2, **89.6%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.

## Your data

- Reads: 49,885 pairs analysed (up to 237 + 232 nt).
- Amplicons: median 426 nt; 99% are 401-427 nt (all: 247-451), measured from the 99% of pairs that merge at full read length.
- Quality scores take only 4 distinct values, as with binned scores on newer Illumina instruments. DADA2's error model can need adjusting for such data.
- Overlap at the recommended settings: about 17 nt on average.

![summary](summary.png)

## How these settings were chosen

- 49,885 of 2,743,627 read pairs were analysed. Totals below 439 nt (your longest amplicon plus the minimum overlap of 12) were not tested.
- Totals were tested from 439 nt in steps of 2 nt, stopping once the share of usable pairs stopped improving. For each, a spread of forward/reverse splits among those with the highest filter pass rate were quality-filtered (maxEE 2) and merged with vsearch. Settings within 1 percentage point of the best yield count as tied; the shortest overlap won.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, and chimeras or pairs that merge at the wrong offset count as usable here, so the yield can be optimistic. Run DADA2 once with these settings to confirm.
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.4.0_
