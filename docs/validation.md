# Validation: the evidence behind each recommendation

AmpliTrim recommends settings, so this page goes through what it recommends and says, for each, what was tested, what held and what did not. The numbers come from real DADA2 runs on 12 public datasets; the scripts, the criteria written down before each run and the full write-up are in this repository (`validation/`, `TEST_REGISTER.md` and `docs/FINDINGS.md`). This page is about the released package, which is in the [software repository](https://github.com/paytonyau/amplitrim). [How it was tested](#how-it-was-tested) is at the end.

## Truncation lengths (`truncLenF`, `truncLenR`)

The total, `truncLenF + truncLenR`, is the shortest that can merge your longest amplicons: the amplicon plus the minimum overlap of 12 nt. How it is split between the two reads, and what it keeps, is searched on your reads.

| claim | result |
|---|---|
| Below that total, merging collapses | Yes. On mock community A, 54 to 55% of pairs merged 10 nt below the floor, 42 to 43% at 20 nt below, 25% at 40 nt below, against 74 to 76% at the floor. On V4-V6 reads no pair merged at three totals below it. Soil V4 is the exception: 10 nt below the floor a cut with a short reverse read merged 67.1%, as many as the best cut at the floor (67.3%); the floor is set to cover all but the longest 0.2% of the amplicons. |
| At that total, the recommendation is close to the best that covers the amplicons | In the development rounds (earlier builds), on five datasets where DADA2 was run on a grid of truncations in steps of 10 nt, it was 0.0 to 0.3 points below the best setting. **Soil V4 is the exception.** With the recommendation of an earlier build (120/118) it was 6.8 points below the best cell (130/108, 59.6% against 66.4%); with the final build's (124/115) it is **4.7 points below** (62.6% against 67.3% for 134/105; criterion of 1 point not met). In both runs the best cells keep the same total, give the forward read more of it and the reverse read less. The soil V4 pick also moved by 11 nt across disjoint sets of samples. |
| It merges as many pairs as FIGARO's truncation, at the same `maxEE` | Equivalent, not better. AmpliTrim uses all read pairs; FIGARO, which needs reads of one length, could use between 8% and 98% of them. |
| The share of pairs that pass DADA2's quality filter is predicted | Within 0.8 points of DADA2's on all 12 datasets (mean absolute error 0.2 points). FIGARO's own predicted retention was off by 0.2 to 7.7 points. |
| The yield estimate ranks settings in DADA2's order | On 9 of the 11 datasets where it was ranked; rank correlation 0.90 and 0.80 on the other two. |

- **On soil V4 the reverse read is what matters. Over-merging, the quality filter and vsearch's view of the overlap do not explain the gap.** In the grid, DADA2's merged share fell as the reverse read got longer: 66.4%, 62.3% and 58.2% for reverse reads of 105, 115 and 125 nt (each the mean of three forward lengths). It hardly changed with the forward read: 61.5%, 62.6% and 62.9% for 114, 124 and 134 nt. At every cut near the floor, 99.9% of the pairs pass the filter. At both 120/118 and 130/108, vsearch merges 99.3% of the reference pairs correctly. So what separates the cells happens after the filter, in denoising or in DADA2's own merging, which AmpliTrim does not model. That is an inference from what was ruled out, not a measurement. AmpliTrim prefers the balanced split among near-ties. A pre-registered test asked whether a forward-heavier split of the same total is better. It was better in 3 of 6 datasets (soil V4 by 4.7 points, mock A by 0.5, V1-V2 by 0.1). This is short of the 4 of 6 that we had set in advance, so we kept the balanced split. The grid used 3,000 pairs per sample; whether the gap persists at full depth was not tested.
- **The overlap is 12 nt, not a larger margin.** That is DADA2's default for `mergePairs`, and real merges are tested at it, where DADA2's author gives a rule of thumb of the amplicon plus 25 nt. No direct comparison with the 25 nt margin was run. In simulation, merges at the wrong length were 0.00% at 12 nt unless a long perfect repeat lay across the overlap; then a longer overlap helped. On real soil V4 reads, against a reference of the same pairs merged from their full-length reads (48,227 pairs), merges at the wrong length were 0.00% at all 18 cuts tested. Extra overlap is not free. On soil V4, each extra 10 nt of overlap added about 0.023 expected errors per pair in the overlap (2.3 per 100 pairs). When the 10 nt went to the reverse read, DADA2's merged share fell by 4.4 to 4.5 points. When the same 10 nt went to the forward read, the share rose by 1.1 and then 0.3 points. Three cuts with the same total, and so the same overlap, merged 58.1%, 62.6% and 67.3%.
- **The yield estimate is optimistic.** It stops before DADA2 denoises the reads. It overshoots the merged share by -2.2 to +6.7 points on the six mock communities and by +3.4 to +7.0 points on four real datasets. On ITS (at 5,000 pairs per sample) and on soil V4, it overshoots by +21.1 and +32.1 points. Use it to compare settings, not as the share of reads you will keep.
- **Merging more strictly does not fix it.** The merge is tolerant of differences in the overlap (up to 10, vsearch's default) because DADA2 merges denoised reads. A strict merge of the raw reads, with no difference allowed, was farther from DADA2: 12.5 points on average, against 6.5 for the tolerant merge. It also ranked the settings worse (rank correlation +0.76 against +0.92). The reason is that, on raw reads, it refuses pairs that denoising would repair. On the six mock communities it undershot by 14.8 points on average. No number of allowed differences from 1 to 7 matched DADA2's merging better than the default of 10 (mean absolute difference 6.6 points at 10 and 27.6 at 0). The part of the overshoot that is left is the loss in denoising, which vsearch cannot see: at least 18.8% of the filtered pairs on ITS and 13.3% on soil V4, against about 1% on the mock communities.

## The coverage setting (the longest 0.2% of the amplicons)

AmpliTrim sets aside the longest 0.2% of the amplicons when it computes the total length, so that a few very long amplicons do not force long, noisy reads on everything. Leaving out more would allow shorter cuts. DADA2's merged share of the input pairs (`maxEE` 2) with the longest 0.2% (the default), 1% and 5% left out:

| dataset | 0.2% | 1% | 5% |
|---|---|---|---|
| D6300 V4-V5 | 57.2% | 89.8% | 89.8% |
| soil V4 | 66.8% | 67.2% | 72.3% |
| V1-V2 | 84.2% | 84.2% | 84.2% |
| inter-laboratory V1-V2 | 95.9% | 95.9% | 93.8% |
| mock community B | 19.9% | 19.9% | 19.9% |
| wastewater V3-V4 | 43.6% | 43.6% | 41.4% |

- **On most datasets, leaving out more changes little or costs.** 1% changed nothing on four of the six, and gained 0.4 points on soil V4; 5% cost about 2 points on two (inter-laboratory V1-V2 -2.1, wastewater -2.2).
- **On two datasets it matters a great deal.** On D6300 V4-V5 about 1% of the amplicons are 550 nt against a median of 374 (a tail suspected to be eukaryotic 18S amplified by the same primers; not checked). The default forced nearly full-length reads (282/280) and merged 57.2%; leaving out 1% gave 268/118 and 89.8%, with 99.9% of the reads in exact ASVs. On soil V4, leaving out 5% gained about 5 points, on 36 of 36 samples.
- **The price is the taxa with the longest amplicons, lost for good.** On mock community B the longest amplicon is a real strain, so leaving out more would lose it. On mock A, leaving out the longest 0.2% instead of 1% recovered one more strain, for 0.2 points of merged pairs.
- The default was not changed, and the command line has no option for it. If your reads contain a tail of off-target amplicons, remove them before running AmpliTrim ([user guide](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#choosing-the-cut-points)).

## `maxEE`

Per read, the limit that covers three quarters of the reads at the recommended cut, never below 2 and never above 5. **This is the least tested recommendation.** The percentile and the cap were chosen after the validation data had been seen, and the rule has not been tested on independent data. In a pre-registered test it did not meet its adoption criteria: it merged about as many pairs as relaxing to `2,5` on average and fewer on three of five datasets. We keep it as a trade-off that we state openly. DADA2's merged share of the input pairs on the five datasets where it differs from the default:

| dataset | `maxEE` 2 | `2,5` | the rule | the rule's `maxEE` |
|---|---|---|---|---|
| mock A | 62.7% | 76.1% | 69.5% | 2,3 |
| mock B | 19.9% | 49.7% | 55.0% | 4,5 |
| inter-laboratory V4-V6 | 41.6% | 60.4% | 57.9% | 3 |
| D6300 V4-V5 | 57.2% | 68.8% | 63.3% | 2,3 |
| wastewater | 43.6% | 51.1% | 57.9% | 4,3 |
| mean of the five | 45.0% | 61.2% | 60.7% | |

It merges more pairs than the default on all five (6.1 to 35.2 points). It merged more than `2,5` on mock B and wastewater, and fewer on mock A (6.7 points), V4-V5 (5.4) and V4-V6 (2.6). On the six clean datasets and ITS it leaves `maxEE` at 2. The cap costs pairs where the reads are very noisy: without it the rule would have gone to 5,9 on mock B and merged 75.3%, as many as FIGARO's pick (74.3%), against 55.0%. On mock B, the one mock where accuracy differed, 73.8% of the reads fell in identical ASVs against 74.7% with `2,5` (the default: 68.7%); on the other two mocks it was as accurate. The rule relaxes each read according to its own noise, so the forward read can be looser than the reverse one: on wastewater the forward limit is 4 and the reverse limit 3. The limit `2,5` relaxes only the reverse read. A percentile of 0.80 matched FIGARO's choices more closely and merged more pairs (68.3% on average), but it lets looser limits through (up to 9 on mock B). We chose 0.75 because it stays closer to DADA2's default.

## ITS (`--its`)

ITS was tested on one dataset (37 spring-wheat samples, PRJEB58189), so it is a case study, not a general result. The recommendation is no fixed truncation length and `maxEE` 2.

- **The filter is replicated exactly.** On the same 49,720 pairs AmpliTrim and DADA2's `filterAndTrim` gave 74.30% of the pairs passing, and 87.47% of the forward reads passing alone. The "where the pairs go" rows and the forward-only share of the report are DADA2's own filter.
- **The merging rows are an upper bound.** AmpliTrim's estimate put 1.9% of the pairs at "pass the filter but do not merge" and the yield at 72%; DADA2 merged 47% of the input pairs. DADA2 lost 15.9 points in denoising, which AmpliTrim does not model, and 11.5 points in merging: 5.1 for an overlap of fewer than 12 matching bases and 6.4 for mismatches in the overlap. No pair was refused for an overhang.
- **`maxEE` is not where ITS loses its pairs.** Relaxing it to `2,5` raised the share passing from 74% to 86% but merged only 2.7 points more.

## What the report's estimates are worth

- **The share passing the filter is exact** (above). **The yield is an upper bound** (above).
- **Accuracy is not assessed.** It was measured on mock communities only, and accuracy in rich real samples was not measured.

## What DADA2's author has said about `maxEE`

Benjamin Callahan, in the repository's issues; each point is from the thread named:

- `maxEE` is the recommended primary quality filter; `truncQ` and `minQ` are not recommended in typical workflows, and the default `truncQ = 2` is permissive on purpose, so it should come into play little (issue #1909, March 2024).
- "maxEE=2 is a good general purpose value. Relaxing it somewhat is also OK", because DADA2's error model uses the quality scores and is "reasonably robust to some higher-noise reads" (issue #2136, July 2025, about full-length PacBio reads).
- The number of reads is not the goal. Keeping about 40% after filtering is fine (said about full-length PacBio reads), and "pushing through more reads of lower quality may not improve the relative abundance estimates, it could even make it worse" (issue #2136). Stricter filtering is fine "when you are just giving up a marginal (e.g. 10%) of reads" (issue #236, November 2020).
- On Illumina reads, trimming is usually more effective than raising `maxEE`, because quality often crashes at the end of the read, especially on 2x300 kits (issue #236, September 2018).
- His own procedure, which he gives "in general", has three steps. First, choose `truncLen` to avoid the worst of the quality profile while keeping enough sequence to overlap healthily. He writes this as forward + reverse `truncLen` above the amplicon length + 25, erring "a bit on the more sequence side". Second, leave `truncQ` at 2. Third, try a couple of `maxEE` values until a satisfactory number of reads passes. This is from issue #236 (April 2017), for a 300 nt amplicon on 2x300 reads, where the cuts he suggested overlapped by about 100 nt. It is a rule of thumb for choosing cuts by eye, not a DADA2 requirement: the default minimum overlap of `mergePairs` is 12.
- He gives no numerical ceiling for `maxEE` in these threads. The numbers sometimes quoted there (2 to 6, or "up to 5 or even 8" for PacBio) were written by the people who opened the issues.

What this means here. Shortening `truncLen` is not an option when AmpliTrim relaxes `maxEE`, because its cut is already the shortest that can merge your longest amplicon plus the minimum overlap. That minimum is 12 nt (DADA2's default for `mergePairs`) and AmpliTrim tests real merges at it, instead of adding a margin chosen by eye as his rule of thumb does (25 nt). And because the author argues that more reads is not the aim, relaxing `maxEE` is worth doing only when the reads need it, which is why the rule leaves clean reads at 2.

## How it was tested

- **Real DADA2 runs on 12 public datasets.** Six are mock communities: V1-V3 twice, V1-V2 and V4-V6 from an inter-laboratory study, and two D6300 mocks at V3-V4 and V4-V5. Six are real datasets: V1-V2, soil V4, wastewater, a spike-in study and swabs at V3-V4, and ITS.
- **Criteria before the run.** Each claim was written down with its margin and verdict words before the run, and the outcome is reported whether it was met or not (`TEST_REGISTER.md` in the archive).
- **Which software.** The runs used the development build that the released package grew from. The released package has the same quality filter and the same search, without the FIGARO comparison and the exploration options, and with the `maxEE` rule built in. DADA2's result depends only on the settings it is given. So we compare the settings that the released package recommends on the 12 datasets with the settings that were run through DADA2, dataset by dataset. This is the release check (section 21 of the register), and its outcome is recorded there.
- **No reads in this repository.** The reads are public, and the accessions are in Table 1 of `docs/FINDINGS.md`. The README of the software shows how to make a small simulated dataset for trying the tool.
- **Reproducing it** needs R with DADA2, vsearch, cutadapt and FIGARO, and the public reads. The steps are in `validation/README.md`. The tests of the released software are in the software repository.
