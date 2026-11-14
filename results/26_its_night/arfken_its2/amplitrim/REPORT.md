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

- About **37.6%** of read pairs should pass quality filtering and merge: 37.6% pass the filter, and 99.9% of those merge.
- maxEE is 2, DADA2's default. AmpliTrim does not adjust it for ITS: reads are not cut at a common position, so there is no cut to measure them at (each read is only trimmed where its quality first falls to 2 or below), and the rule has not been tested on reads of different lengths.

## Things to check

- Only 78% of read pairs could be merged at full read length (amplicons longer than about 588 nt, or some reads shorter than the rest); the lengths above come from the pairs that merged.
- Only 38% of read pairs are expected to pass filtering and merge.

## Where the pairs go

Of the 36,721 read pairs analysed, with maxEE 2 and a minimum length of 50 nt (the filter rows replicate DADA2's filter; the merging rows are vsearch's estimate):

| step | share of the pairs |
|---|---|
| too short after the cut at the first base of quality 2 or below (under 50 nt) | 0.0% |
| contain an N | 0.0% |
| more expected errors than maxEE 2 | 62.4% |
| **pass the filter** | **37.6%** |
| pass the filter but do not merge | 0.0% |
| **pass the filter and merge** | **37.6%** |

The first rows are DADA2's own filter, replicated. The last two are vsearch's estimate on the raw reads, which cannot see DADA2's denoising or its own merging (see Limits).
- **Using the forward reads alone,** 76.5% of the reads would pass the same filter, against 37.6% of the pairs that pass and merge. DADA2's author says that using only the forward reads is "totally reasonable", and that when the reverse reads are bad enough they "will cost more in sensitivity to low-frequency stuff than they add in a lower FP rate" (DADA2 issue #236, 2017). AmpliTrim does not measure what this does to your ASVs.

## maxEE: retention against read quality

ITS reads are not given a fixed truncation length, so maxEE is the one setting that trades reads for quality: each read is still cut at its first base of quality 2 or below, but a read with a noisy tail is not shortened, it is dropped when its expected errors exceed maxEE. For each value: the share of pairs that pass DADA2's filter, the estimated yield, how long the kept reads are, and how noisy they are.

| maxEE | pairs passing the filter | estimated yield | mean length kept (F / R) | expected errors per kept read (F / R) | expected errors in the reads added beyond maxEE 2 (F / R) |
|---|---|---|---|---|---|
| 1 | 17.5% | 17.5% | 278 / 278 | 0.16 / 0.60 | - |
| **2 (used)** | 37.6% | 37.6% | 277 / 278 | 0.25 / 1.06 | - |
| 2,3 | 51.1% | 51.0% | 277 / 278 | 0.32 / 1.43 | 0.50 / 2.47 |
| 2,4 | 60.1% | 59.9% | 277 / 278 | 0.37 / 1.74 | 0.58 / 2.87 |
| 3 | 51.8% | 51.8% | 277 / 278 | 0.35 / 1.44 | 0.61 / 2.45 |
| 2,5 | 66.3% | 65.5% | 277 / 278 | 0.42 / 1.99 | 0.64 / 3.21 |
| 3,5 | 68.8% | 67.8% | 277 / 278 | 0.49 / 2.05 | 0.78 / 3.23 |
| 5 | 70.1% | 68.8% | 277 / 278 | 0.55 / 2.08 | 0.90 / 3.25 |

At maxEE 2, **37.6%** of the read pairs pass DADA2's quality filter at these truncation lengths.
DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK. DADA2's error model uses the quality scores, so it tolerates some noisier reads. But pushing through more lower-quality reads may not improve the relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: the expected errors in the kept reads.
A low share is not a problem in itself: DADA2's author says that the number of reads is not the aim, and that keeping about 40% after filtering is fine (DADA2 issue #2136, about PacBio reads).

The yield is an upper bound: it stops before DADA2's denoising and chimera removal, so DADA2 keeps less, and the less so at low depth or high diversity (on a spring-wheat ITS study DADA2 merged 0.65 of the estimate at 5,000 pairs per sample and 0.87 at full depth). The gain from a looser maxEE is real but smaller than the estimate promises, and AmpliTrim does not measure what a looser maxEE does to the accuracy of your results.

## Your data

- Reads: 36,721 pairs analysed (up to 300 + 300 nt).
- Amplicons: median 318 nt; 99% are 266-399 nt (all: 110-419), measured from the 78% of pairs that merge at full read length.
- Overlap at the recommended settings: about 209 nt on average.

![summary](summary.png)

## How these settings were chosen

- 36,721 of 73,441 read pairs were analysed.
- No fixed truncation length is tested, as in DADA2's ITS workflow: each read is cut at its first base of quality 2 or below (`truncQ = 2`), then the pair is kept if both reads are at least 50 nt, free of Ns and within maxEE 2. The kept pairs are merged with vsearch.

## Limits

- The estimates come from merging raw reads with vsearch, not from DADA2. Run DADA2 once with these settings to confirm.
- Pairs whose amplicon is longer than the reads can span cannot merge, with any setting.
- The 'pass the filter but do not merge' row can be far too low: DADA2 also loses pairs in denoising and in its own merging, and an estimate from the raw reads sees neither (the user guide gives the figures for the one ITS dataset we checked).
- Only read quality and mergeability are assessed, not taxonomic accuracy.
- The yield is measured on the raw reads, before DADA2 denoises them. In diverse samples (soil, for example) at modest depth DADA2 can keep far fewer reads than this estimate. Use it to compare settings, not as the final retained fraction.
- The reads are assumed to have adapters and primers removed and no quality trimming.

_AmpliTrim 0.2.0_
