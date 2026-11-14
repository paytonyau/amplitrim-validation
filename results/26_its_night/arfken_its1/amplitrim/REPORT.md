# AmpliTrim report

## Recommended settings

Do not truncate. ITS amplicons differ in length, so a fixed `truncLen` would discard every read shorter than it.

| truncLenF | truncLenR | maxEE | minLen | minOverlap |
|---|---|---|---|---|
| **0 (none)** | **0 (none)** | 2 | 50 | 12 |

In DADA2 (R):

```r
out <- filterAndTrim(fwd, filtF, rev, filtR, maxN = 0,
                     maxEE = c(2, 2), truncQ = 2, minLen = 50)
# ... learnErrors(), dada() ...
merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = 12, maxMismatch = 0)
```

In QIIME 2:

```
qiime dada2 denoise-paired ... --p-trunc-len-f 0 --p-trunc-len-r 0 \
    --p-max-ee-f 2 --p-max-ee-r 2
```

## What to expect

- About **61.9%** of read pairs should pass quality filtering and merge: 62.0% pass the filter, and 99.9% of those merge.
- maxEE is 2, DADA2's default. AmpliTrim does not adjust it for ITS: reads are not cut at a common position, so there is no cut to measure them at (each read is only trimmed where its quality first falls to 2 or below), and the rule has not been tested on reads of different lengths.

## Things to check

- Only 84% of read pairs could be merged at full read length (amplicons longer than about 588 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 62% of read pairs are expected to pass filtering and merge.

## Where the pairs go

Of the 37,076 read pairs analysed, with maxEE 2 and a minimum length of 50 nt (the filter rows replicate DADA2's filter; the merging rows are vsearch's estimate):

| step | share of the pairs |
|---|---|
| too short after the cut at the first base of quality 2 or below (under 50 nt) | 0.0% |
| contain an N | 0.0% |
| more expected errors than maxEE 2 | 38.0% |
| **pass the filter** | **62.0%** |
| pass the filter but do not merge | 0.1% |
| **pass the filter and merge** | **61.9%** |

The first rows are DADA2's own filter, replicated. The last two are vsearch's estimate on the raw reads, which cannot see DADA2's denoising or its own merging (see Limits).
- **Using the forward reads alone,** 83.2% of the reads would pass the same filter, against 61.9% of the pairs that pass and merge. DADA2's author says that using only the forward reads is "totally reasonable", and that when the reverse reads are bad enough they "will cost more in sensitivity to low-frequency stuff than they add in a lower FP rate" (DADA2 issue #236, 2017). AmpliTrim does not measure what this does to your ASVs.

## maxEE: retention against read quality

ITS reads are not given a fixed truncation length, so maxEE is the one setting that trades reads for quality: each read is still cut at its first base of quality 2 or below, but a read with a noisy tail is not shortened, it is dropped when its expected errors exceed maxEE. For each value: the share of pairs that pass DADA2's filter, the estimated yield, how long the kept reads are, and how noisy they are.

| maxEE | pairs passing the filter | estimated yield | mean length kept (F / R) | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|---|
| 1 | 51.4% | 51.4% | 196 / 196 | 0.10 / 0.19 | - |
| **2 (used)** | 62.0% | 61.9% | 204 / 205 | 0.17 / 0.40 | - |
| 2,3 | 68.4% | 68.2% | 209 / 210 | 0.21 / 0.59 | 0.61 / 2.47 |
| 2,4 | 72.9% | 72.6% | 213 / 214 | 0.24 / 0.77 | 0.64 / 2.88 |
| 3 | 69.2% | 69.0% | 210 / 210 | 0.24 / 0.61 | 0.80 / 2.42 |
| 2,5 | 76.0% | 75.4% | 215 / 216 | 0.27 / 0.92 | 0.68 / 3.23 |
| 3,5 | 77.9% | 76.9% | 216 / 217 | 0.32 / 0.98 | 0.89 / 3.22 |
| 5 | 78.9% | 77.5% | 216 / 217 | 0.36 / 1.01 | 1.05 / 3.23 |

At maxEE 2, **62.0%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.
A low share is not a problem in itself: DADA2's author says that the number of reads is not the aim, and that keeping about 40% after filtering is fine (DADA2 issue #2136, about PacBio reads).

The yield is an upper bound: it stops before DADA2's denoising and chimera removal, so DADA2 keeps less, and the less so at low depth or high diversity (on a spring-wheat ITS study DADA2 merged 0.65 of the estimate at 5,000 pairs per sample and 0.87 at full depth). The gain from a looser maxEE is real but smaller than the estimate promises, and AmpliTrim does not measure what a looser maxEE does to the accuracy of your results.

## Your data

- Reads: 37,076 pairs analysed (up to 300 + 300 nt).
- Amplicons: median 215 nt; 99% are 154-395 nt (all: 25-422), measured from the 84% of pairs that merge at full read length.
- Overlap at the recommended settings: about 192 nt on average.

![summary](summary.png)

## How these settings were chosen

- 37,076 of 74,151 read pairs were analysed.
- No fixed truncation length is tested, as in DADA2's ITS workflow: each read is cut at its first base of quality 2 or below (`truncQ = 2`), then the pair is kept if both reads are at least 50 nt, free of Ns and within maxEE 2. The kept pairs are merged with vsearch.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. Run DADA2 once with these settings to confirm.
- Pairs whose amplicon is longer than the reads can span cannot merge, with any setting.
- The 'pass the filter but do not merge' row can be far too low: DADA2 also loses pairs in denoising and in its own merging, and an estimate from the raw reads sees neither (the user guide gives the figures for the one ITS dataset we checked).
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
