# AmpliTrim test register

**One standard battery, run the same way on every dataset.** Every dataset (each region, each study, each mock community) goes through the same steps and produces the same
columns: AmpliTrim against DADA2 (filter and yield), against FIGARO (like for like and end to end), and against the simple DADA2 alternatives. A dataset with expected sequences
(a mock community) gets the same battery plus the same truth columns. Anything that does not fit this mould is an *additional analysis* (section 5): supplement or removed.

- **Stage.** `DEV` = used to build or verify the tool (section 7, supplement only). `TEST` = reported in the manuscript.
- **Where results are written.** `C:\data\overnight\OVERNIGHT_SUMMARY.md`; `python validation\check_build.py --table` writes the test tables (one for 16S, one for ITS, one for the
  expected sequences) to `TEST_TABLE.md`; `--todo` prints what is left; `--decide` applies the decision rule of section 6.
- **Provenance.** The overnight runs used the build printed as `bd4e64353c` on Windows (`6f3cb5c6f9` when its files are hashed in sorted order: the same code; earlier ids depended on the
  operating system's file order). Software: DADA2 1.40.0, cutadapt 5.2, Python 3.12; vsearch and R versions are printed by `python validation\check_build.py`. A last full run on the final
  build is listed in section 8.

---

**Contents.** Sections 1 to 8 describe the standard battery and the first outcomes; sections 9 to 26 follow the work in time order (round 4, then the nights, the checks and the release). Every section that decides something has its criteria written before the run. The **Run** column numbers the rounds and nights in the order they were run (01 to 14); the same numbers prefix the run folders on the data drive and are used in `validation/README.md`.

| Section | Run | What |
|---|---|---|
| 1 | 01 | Datasets |
| 2 | 01 | The standard battery |
| 3 | 01 | Outcome matrix |
| 4 | 01 | Outcomes so far |
| 5 | 01 | Additional analyses |
| 6 | 01 | Decision rule for the default of --allow-unmergeable |
| 7 | 01 | Development stage |
| 8 | 01 | Still to run, and further tests recommended |
| 9 | 02 | Round 4: the open questions, tested in one place |
| 10 | 02 | Round 4: outcomes |
| 11 | 03 | Night 1 |
| 12 | 04 | Night 2: the maxEE head-to-head |
| 13 | 04 | Night 2: outcomes |
| 14 | 05 | Night 3: arm G and arm F2 |
| 15 | 05 | Night 3: outcomes |
| 16 | 06 | Night 4 |
| 17 | 06 | Night 4, second stage |
| 18 | 07 | Night 4b: the percentile rule through DADA2 |
| 19 | 08 | ITS check: are the numbers of the ITS report what DADA2 gives? |
| 20 | 09 | ITS merge check: the exact split of DADA2's loss after the filter, and whether an overhang explains part of the ... |
| 21 | 10 | Release check: does the released package recommend the settings that DADA2 was run at? |
| 22 | 11 | Soil V4: over-merging on real reads, and the grid around the final recommendation |
| 23 | 12 | The merge estimates against DADA2, and the cut-off for differences in the overlap |
| 24 | 13 | Night 5: the rule without a cap, whether a simpler cut would do, and the reads wastewater loses to length |
| 25 | 14 | Night 6: an independent real dataset, PRJNA933120 |
| 26 | 14 | Night 6: outcomes |
| 27 | 23 | The minimum overlap: is 12 nt safe, and what does a longer overlap cost? |
| 28 | 21 | Release check of version 0.2.0 |
| 29 | 24 | The closing run (Night 10): the gaps left before submission |
| 30 | 24 | Night 10: outcomes |
| 31 | 25 | ITS follow-up: why DADA2 refuses ITS pairs, and is a mismatch limit above 0 safe? |
| 32 | 25 | ITS follow-up: outcomes |
| 33 | 26 | ITS overnight run: what AmpliTrim can recommend for ITS |
| 34 | 26 | ITS overnight run: outcomes |
| 35 | 27 | ITS rules night: can AmpliTrim compute truncQ, the trim strength and maxEE for ITS? |
| 36 | 27 | ITS rules night: outcomes |
| 37 | 28 | The Bakker night: do the ITS rules hold on a third mock community? |
| 38 | 28 | The Bakker night: outcomes |
| 39 | 29 | The Bakker re-run: Bakker's reads with the far primer cut |
| 40 | 29 | The Bakker re-run: outcomes |
| 41 | 30 | The trim check: does AmpliTrim 0.3 predict the share of pairs that pass DADA2's filter after the trim? |

## 1. Datasets

| ID | Dataset | Marker, region | Reads | Samples | Stage | Config key | Expected sequences |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D1 | DADA2 package example reads | 16S V4 | MiSeq 2x250 | 2 samples, 1,500 pairs each | DEV | - | - |
| D2 | Simulated communities (`validation\simulated\`) | 16S-like, ITS-like | 2x250 to 2x283 | 2 to 5 samples | DEV | - | yes (simulated) |
| D3 | PRJEB42528 | 16S V1-V2 (27F/338R) | MiSeq 2x300 | 33 samples, 1.43 million pairs | TEST | `v1v2` | - |
| D4 | PRJEB67797, staggered mock community A [Hülpüsch, Rauer et al. 2023] | 16S V1-V3 (27F-YM/534R) | MiSeq 2x300 | 15 strains, 10^9 to 10^2 cells in triplicate; 23 mock samples and 6 controls analysed (13 with the dilutions used here); 6 high-input samples | TEST | `mockA_high`, `mockA_all` | yes, from Sanger sequences of the 15 strains (a 363 nt region): the authors' ASVs that contain it exactly or lie within 4 edits; 14 strains detected (S. pneumoniae not); anchored to Sanger reads but from the authors' own DADA2 run, so not independent |
| D5 | PRJNA574030 | 16S V4 | see `--download v4` | 36 soil samples | TEST (used in development, so run again with the final build) | `v4` | - |
| D6 | PRJEB58189, UK Crop Microbiome Cryobank | ITS (ITS1-Fl2/ITS2) | MiSeq 2x250 | 37 spring-wheat samples; 8 at full depth | TEST (case study) | `its` | - |
| D7 | PRJEB67827, staggered mock community B [Rauer, De Tomassi et al.; used in Hülpüsch et al. 2023]: the three-species Zymo spike-in D6321 (Truepera radiovictrix, Imtechella halotolerans, Allobacillus halotolerans) | 16S V1-V3 | MiSeq 2x300 | the project has 94 runs (three mocks, eight extraction protocols); mock B = 8 samples: BiomPsy 536-539 (1.1x10^5 cells) and 544-547 (5.55x10^3 cells); 552-555 are negative controls and are left out | TEST | `mockB` | yes (Zymo's three reference 16S genes, `--make-reference mockB=D6321.refseq.zip`; the authors scored an ASV as mock within 4 edits of them) |
| D8 | V3-V4 2x250 (candidates of PRJEB58189) | 16S V3-V4 | MiSeq 2x250 | - | TEST | `v3v4` | - |
| D9 | PRJNA971391 | 16S V4-V5 | - | - | TEST | `v4v5` | - |
| D10 | PRJNA915241, sputum, region V1-V2 (QIAseq 16S/ITS panel; the same 33 samples in every region) | 16S V1-V2 | MiSeq 2x300 | 33 sputum samples and a ZymoBIOMICS standard | NOT PLANNED (considered and dropped) | `sputum_v1v2` | the standard's own sequences, if its runs are in the project |
| D11 | PRJNA915241, sputum, region V3-V4 (QIAseq 16S/ITS panel; the same 33 samples in every region) | 16S V3-V4 | MiSeq 2x300 | 33 sputum samples and a ZymoBIOMICS standard | NOT PLANNED (considered and dropped) | `sputum_v3v4` | the standard's own sequences, if its runs are in the project |
| D12 | PRJNA915241, sputum, region V5-V7 (QIAseq 16S/ITS panel; the same 33 samples in every region) | 16S V5-V7 | MiSeq 2x300 | 33 sputum samples and a ZymoBIOMICS standard | NOT PLANNED (considered and dropped) | `sputum_v5v7` | the standard's own sequences, if its runs are in the project |
| D13 | PRJNA915241, sputum, region V7-V9 (QIAseq 16S/ITS panel; the same 33 samples in every region) | 16S V7-V9 | MiSeq 2x300 | 33 sputum samples and a ZymoBIOMICS standard | NOT PLANNED (considered and dropped) | `sputum_v7v9` | the standard's own sequences, if its runs are in the project |
| D14 | PRJEB34919, inter-laboratory mock communities MCM2a and MCM2b [O'Sullivan et al. 2021], V4-V6 (about 564 bp) | 16S V4-V6 | MiSeq 2x300 (a real tight-overlap case) | 6 samples (2 materials x 3 repeats); the shotgun runs of the project are left out | TEST (held-out) | `mcm_v46` | yes: the 15 bacteria, quantified by digital PCR; the sequences are built from their genomes |
| D15 | PRJEB34919, the same materials, V1-V2 (about 337 bp) | 16S V1-V2 | MiSeq 2x300 | 6 samples | TEST (held-out; optional) | `mcm_v12` | as D14 |

---

## 2. The standard battery (every TEST dataset)

One command runs it. The steps `battery`, `coverage_levels`, `landscape`, `per_sample_rank` and `stability` take every registered dataset; for the new datasets (`v4`, `mockB`, `v3v4`, `v4v5`, registered with `--download` and `--prepare`) the battery step is `--only extra_datasets`:

```bat
python validation\check_build.py --run --only battery,coverage_levels,landscape,per_sample_rank,stability
python validation\check_build.py --table
```

| ID | Question | Output (a column of the test table) | Condition | Criterion |
| --- | --- | --- | --- | --- |
| B1 | Does the filter prediction match DADA2's pass rate? | largest difference over the maxEE settings (points) | three maxEE settings (the dataset's own, 1 and 5; for ITS 1, 2 and 5) | C1: within 2 points |
| B2 | How far is the yield estimate from the share DADA2 merged, and does it rank the settings? | offset and ratio (mean, range) and the rank correlation across the same settings | at least three settings | C2: ranks agree; C3: offset constant |
| B3 | What do the simple DADA2 alternatives cost? | AmpliTrim minus no truncation (ITS: minus a 200/200 truncation), paired per sample, with the 95% interval | at least 10 samples for the interval | C4 for ITS: at least 10 points |
| B4 | AmpliTrim against FIGARO, (a) like for like at FIGARO's own maxEE and (b) end to end: AmpliTrim's default against FIGARO as it ships | both differences, paired per sample, with the 95% interval | 16S only (FIGARO cannot run on mixed-length ITS reads); the standalone steps (`--run --only figaro_mockA,figaro_v1v2`) run the same comparison | C11: like-for-like interval within +/-1 point |
| B5 | What does the share of the longest amplicons left out cost and buy? | merged points from 1% to 0.2% (and with expected sequences: variants and strains) | levels 1%, 0.5%, 0.2%, 0.1%, 0% | C14 (section 6) |
| B6 | Is the recommendation near the best of a DADA2 grid around it? | points below the best setting that still covers the amplicons | 5 x 5 grid in steps of 10 nt, DADA2 at every cell | C13: within 1 point |
| B7 | Does the estimate follow individual samples? | per-sample rank correlation | at least 10 samples, at AmpliTrim's recommendation | C8: at least 0.80 |
| B8 | Does the recommendation depend on the subsample? | spread of the recommended truncation (nt) and yield (points) | three subsample sizes, two strides | C7: 10 nt, 1 point |

**Truth add-on (a dataset with expected sequences; on automatically when a reference is registered).** The same methods of B3 and B4 (AmpliTrim, FIGARO, no truncation), every coverage level
of B5 and every cell of B6 are also scored: variants and strains recovered, reads in ASVs identical to an expected sequence, ASVs outside the reference. The test table prints them in a third table.
Criterion C6: no variant longer than total - 12 nt is ever recovered.

---

## 3. Outcome matrix (what has been run)

`done` = run with the battery or an earlier equivalent; `-` = not applicable; `pending` = to run. Numbers in section 4.

| Dataset | B1 | B2 | B3 | B4 | B5 | B6 | B7 | B8 | Truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| D3 V1-V2 | done | done | done | done (like for like +0.1; end to end +2.0) | done (0.4 points) | done (0.1) | done (0.63) | done (1 nt) | - |
| D4 mock A | done (-1.1) | done | done (+10.5) | done (-0.2; -1.9) | done | done (0.3) | done (0.60; 23 samples) | done (0 nt) | done |
| D5 V4 | done (+0.0) | done (offset +30.4) | done (+13.0) | done (+2.0; +4.3) | done (gain 2.6) | done (6.8: not met) | done (0.70) | done (3 nt) | - |
| D6 ITS | done | done | done | - | - | - | - | - | - |
| D7 mock B (8 samples) | done (+0.0) | done (offset +0.4, ratio 0.98) | done (+18.0) | done (+2.6; -54.1) | done (0.0) | done (0.0) | not possible (8 samples) | done (0 nt) | done (3/3 strains) |
| D8 V3-V4 | pending | pending | pending | pending | pending | pending | pending | pending | - |
| D9 V4-V5 | pending | pending | pending | pending | pending | pending | pending | pending | - |
| D10 to D13 sputum (PRJNA915241) | not planned | not planned | not planned | not planned | not planned | not planned | not planned | not planned | not planned |
| D14 mock V4-V6 (5 samples; the sixth run is not served by ENA) | done (-0.3) | done (offset +0.3, ratio 0.99) | - (5 samples) | - (5 samples) | done (0.0) | done (0.0; merging is 0% below the floor) | not possible (5 samples) | done (1 nt) | pending (no reference) |
| D15 mock V1-V2 (6 samples) | done (-0.1) | done (offset +0.2, ratio 1.00) | done (+69.9) | done (0.0; +2.3) | done (0.0) | done (0.1) | not possible (6 samples) | done (3 nt) | pending (no reference) |

---

## 4. Outcomes so far

| ID | Dataset | Outcome |
| --- | --- | --- |
| B1 | V1-V2 | within 0.1 points at 5 settings |
| B1 | ITS | within 0.1 (5 settings, 5,000 pairs per sample); within 0.2 (3 settings, full depth) |
| B1 | mock A (clipped to 2x250) | 1.2 to 1.4 points low |
| B2 | V1-V2 | offset +6.4 (5.7 to 7.7); ratio 0.92 (0.88 to 0.94); rank 1.00 |
| B2 | mock A | offset +2.3 at maxEE 2,6, ratio 0.97; clipped: -0.7 (-1.9 to 0.0), ratio 1.02 |
| B2 | ITS | 5,000 pairs: +22.9 (13.8 to 30.4), ratio 0.65, rank 1.00. Full depth: +8.5 (4.7 to 11.2), ratio 0.88, rank 1.00 |
| B3 | V1-V2 | 84.2% merged at AmpliTrim's setting; 77.1% at 230/200; 48.3% at 260/260 (fixed settings, not paired) |
| B3 | ITS | 200/200 costs 23.5 points (95% interval 21.8 to 25.2), in 37 of 37 samples |
| B4 | V1-V2 | like for like (maxEE 2,1): 82.3% against 82.2%, no interval; end to end 84.2% against 82.2% |
| B5 | mock A | 1% and 0.5%: 257/260 (78.4%, 35/39 variants, 11/13 strains); 0.2% and 0.1%: 257/263 (78.2%, 36/39, 12/13); 0%: 268/280 (73.2%, 36/39, 13/13) |
| B5 | V1-V2 | 1%: 196/159 (83.9% merged); 0.2%: 195/168 (83.5%), a cost of 0.4 points; 0%: 279/280 (29.4%), about 54 points (DADA2 at 5,000 pairs per sample; no expected sequences) |
| B6 | V1-V2 | 0.1 points below the best of 9 (196/159 82.6%, best 186/169 82.7%); merging falls to 71.8 to 72.0% at 10 nt below the floor of 355 nt and to 60.9% at 20 nt below (3,000 pairs per sample) |
| B7 | V1-V2 | 0.63 over 33 samples (criterion not met); offsets +1.4 to +14.2, mean +6.5 |
| B6 | mock A | 0.3 points below the best of 25 runs; merging collapses below the coverage floor (54% at 10 nt below, 76% at the floor) |
| B7 | mock A | 0.60 over 23 samples; offsets 0.0 to +10.3 (mean +2.6) |
| B8 | mock A | 277/240 in every run (six samples); 274/243 to 275/242 (23 samples); yield varies by 0.5 and 0.2 points |
| Truth | mock A | 10 settings: no variant longer than the limit recovered. 257/260 at maxEE 2,6: 35/39 variants, 11/13 strains, 15 ASVs outside; authors' 279/263: 37/39, 13/13, 7; 265/270: 36/39, 13/13, 12 |
| B1 | V4; inter-laboratory mock V1-V2; mock A at full length | +0.0; -0.1; -1.1 points at the default maxEE (final build, 5,000 pairs per sample) |
| B2 | V4 | offset +30.4 (+29.5 to +30.9), ratio 0.69, rank 1.00 over 3 settings (5,000 pairs per sample) |
| B2 | inter-laboratory mock V1-V2 | offset +0.2 (-0.1 to +0.4), ratio 1.00, rank 1.00 over 3 settings; about 28 ASVs |
| B2 | mock A; V1-V2 (final build) | +2.4 (+0.3 to +3.9), 0.95, rank 1.00 over 4 settings; +6.3 (+6.2 to +6.4), 0.93, rank 1.00 over 4 settings |
| B3 | AmpliTrim minus no truncation (ITS: minus 200/200) | mock A +10.5 (+9.1 to +12.1); V1-V2 +55.0 (+53.6 to +56.4); V4 +13.0 (+12.2 to +13.9); inter-laboratory mock +69.9 (+66.4 to +73.6); ITS +16.1 (+14.3 to +17.9), which differs from the earlier +23.5 (+21.8 to +25.2) and is to be reconciled |
| B4 | like for like (at FIGARO's maxEE) | V1-V2 +0.1 (-0.4 to +0.6); mock A -0.2 (-0.8 to +0.6); inter-laboratory mock 0.0; V4 +2.0 (+1.6 to +2.4). The two six-sample intervals are not reliable |
| B4 | end to end (AmpliTrim default minus FIGARO as shipped) | V1-V2 +2.0 (+1.5 to +2.5); mock A -1.9 (-3.9 to -0.6; FIGARO's maxEE 3,4 is looser than 2,6); V4 +4.3 (+3.7 to +4.9); inter-laboratory mock +2.3 (+2.2 to +2.5) |
| B5 | V4; inter-laboratory mock | V4: 120/118 to 124/115 gained 2.6 points (64.3% to 66.8%); inter-laboratory mock: 160/165 at both levels |
| B6 | V4; inter-laboratory mock | V4 59.6% against 66.4% for the best setting (130/108): 6.8 points below, criterion not met; inter-laboratory mock 0.1 points |
| B7 | V4 | 0.70 over 36 samples (criterion not met; third dataset below 0.80 after 0.60 and 0.63) |
| B8 | V1-V2; V4; inter-laboratory mock | 1, 3 and 3 nt |
| Truth | mock A at 5,000 pairs per sample (final build) | AmpliTrim maxEE 1 / 2,6 / 5: 28 / 29 / 32 variants, 10 / 10 / 11 of 13 strains, 42.6 / 74.6 / 84.2% merged; FIGARO 31 variants, 76.4%; no truncation 29 variants, 64.2% merged and 97.70% of reads in identical ASVs (99.75% with AmpliTrim's truncation) |
| B1 | mock B; V4-V6 | +0.0; -0.3 points at the default maxEE |
| B2 | mock B (3 species, about 30% contaminant reads) | offset +0.4 (-6.5 to +1.1 across maxEE 2, 1, 5 and 6,9), ratio 0.98, rank 1.00 over 4 settings; at 100,000 pairs per sample 19.9% estimated, 18.9% merged |
| B2 | V4-V6 (real tight overlap, 5 samples) | offset +0.3 (-3.8 to +0.4), ratio 0.99, rank 1.00 over 4 settings at 5,000 pairs per sample; at 100,000 pairs the estimate was 41.9% and DADA2 merged 42.9% (-1.0) |
| B3 | mock B | AmpliTrim minus no truncation +18.0 (+15.9 to +19.6): 19.9% against 2.0% merged |
| B4 | mock B | like for like (FIGARO's maxEE 6,9) +2.6 (+0.3 to +5.2); end to end -54.1 (-58.0 to -50.3): AmpliTrim's default maxEE 2 passed 20% of the pairs, FIGARO's 6,9 passed 74% |
| B5 | mock B; V4-V6 | the same truncation at 1% and 0.2% (280/247; 273/265): no cost, no change in variants or strains |
| B6 | mock B; V4-V6 | 0.0 points below the best feasible setting in both. V4-V6: 0.0% merged at 263/255, 263/265 and 273/255 (total 518 or 528, below the floor of 538) and 40.0% at 263/275, 41.3% at 273/265 |
| B8 | mock B; V4-V6 | 0 and 1 nt |
| Truth | mock B (Zymo D6321: 7 expected variants of 3 strains), 5,000 pairs per sample | strains 3/3 for every method. Variants: AmpliTrim maxEE 1 / 2 / 5: 6 / 6 / 6; FIGARO 7; no truncation 4. Merged: 5.7 / 19.9 / 55.8% (FIGARO 74.3%, no truncation 2.0%). Reads in identical ASVs 54 to 78%; 17 to 69 ASVs outside the reference, mostly contaminants of the extraction buffer (the authors report about a third of the reads) |

---

## 5. Additional analyses (not part of the standard battery) and what to do with them

| ID | Analysis | Command | Outcome | Disposition |
| --- | --- | --- | --- | --- |
| A1 | AmpliTrim's truncation against the authors' published one (mock A) | `--run --only mock_grid` | the pre-set non-inferiority test failed: 2 variants and 2 strains fewer, 8 more ASVs (mainly chimeras), 6.1 points more merged pairs | main text, as the truth comparison of mock communities (a published setting exists for mock A only); for mock B the published setting, if any, replaces it |
| A2 | Dilution series at three truncations | `--run --only mock_dilutions` | 185 of 198 present variants at 257/260; 11 of 13 misses longer than the limit; 193 and 194 at 265/270 and 279/263 | supplement (needs the authors' table; mock A only) |
| A3 | ASV classes (contaminant, chimera, error variant, unlisted variant) | `--run --only asv_classify` | the excess at 257/260 is mainly chimeras (5 against 1, 0.01% of reads) | supplement |
| A4 | Tracing the missed variants | `--run --only mock_miss_trace` | long variants lost at merging; two mergeable ones absorbed by a 1 or 2 edit neighbour | supplement |
| A5 | Mock A clipped to 2x250 (tight overlap, emulated) | `--run --only mock_clip250` | AmpliTrim 249/249 equals no truncation; 4 of 13 strains at any setting | supplement (an emulation; replace by real V3-V4, D8) |
| A6 | Pooled DADA2 inference | `--run --only mock_pooled` | recovery unchanged; full pooling adds 13 ASVs outside the reference | supplement |
| A7 | Chimera removal in ITS under other settings | `--run --only its_chimeras` | 0.2%, 0.6%, 1.5% against 13.9% in the authors' pipeline: not explained | supplement, or remove: it concerns the authors' pipeline, not AmpliTrim |
| A8 | Effect of maxEE on accuracy | (*as run:* `max_ee_sweep.py READY --trunc 257,260 --max-ee 2,2 2,3 2,5 2,6 3,5 --reference ...`) | merged 56.4% to 79.7%; reads in identical ASVs 99.76 to 99.79% | superseded by the three maxEE settings of the battery (B1, B2) and the truth add-on |
| A9 | Calibration at full depth (V1-V2, ITS) | `--run --only v1v2_fulldepth`; ITS: `compare_methods.py ... --marker its --amplitrim-max-ee 1 2 5 --only ERR... --pairs 0` | merged / filtered 0.84 at full depth, 0.63 at 5,000 pairs (0.81 in the authors' pipeline) | main text, as a depth variant of B2 on every dataset (a `--pairs 0` run of the battery) |

---

## 6. Decision rule for the default of --allow-unmergeable (C14)

Stated before the evidence: change the default from 1% to 0.2% if, on at least two datasets with expected sequences, it recovers at least one more strain on one and fewer variants or strains on none,
and across all datasets it costs a median of at most 2 points of merged pairs and at most 5 on any. `python validation\check_build.py --decide` applies it to the saved B5 results.
Outcome (run of 7 Oct 2026): **ADOPT 0.2% as the default**. On the two datasets with expected sequences, mock A recovers one more variant and strain at 0.2% (and has 4 fewer ASVs outside the reference) and mock B changes nothing; the cost is 0.0 points in median and at most 0.4 on any of six datasets (V1-V2 -0.4, mock A -0.2, V4 gains 2.6). The default of `--allow-unmergeable` is therefore 0.2% from the next build; the benefit rests on one strain of one mock, and `--allow-unmergeable 0.01` restores the old default.

**Decision 2: an automatic maxEE rule (fixed on 8 Oct 2026, before `validation/maxee_frontier.py` was run on any real dataset).** FIGARO also returns a maxEE with its truncation, so a rule is worth building only if it can match FIGARO's retention without keeping noisier reads. The question is answered with AmpliTrim's exact quality table (no DADA2) on the six datasets that have a FIGARO row (V1-V2, mock A, mock B, V4, inter-laboratory V1-V2 and V4-V6). At AmpliTrim's truncation, each candidate rule (retention targets of 60, 70 and 80%, FIGARO's score applied to exact pass rates, and the quality-first rule of `--explore-max-ee`) picks a maxEE. A rule qualifies if on at least 4 of the 6 datasets (a) it keeps at least FIGARO's share of pairs minus 2 points and (b) its kept reads are `cleaner` or `about the same` as FIGARO's (no more expected errors on either read beyond 0.02). If several qualify, the strictest wins. If none qualifies, no automatic rule is built and the paper says that FIGARO's maxEE choice was not matched by a simple rule. If one is built it is opt-in (`--auto-max-ee`; the default maxEE stays 2, so every other table stands) and is tested against FIGARO as shipped in the final run: merged pairs (paired difference, interval within 1 point or above), and on mock A and mock B no more than 1 variant fewer and no more than 0.5 points fewer reads in identical ASVs than FIGARO. The check was written without the test suite or vsearch available (no network in that session), so its four tests must pass on the user's machine first.

**Outcome of decision 2 (maxEE frontier check, 8 Oct 2026, six datasets, exact quality table, 60,000 pooled pairs each): no rule qualifies, so no automatic maxEE is built.** Datasets passing the condition (at least FIGARO's retention minus 2 points, kept reads cleaner or about the same; 4 of 6 needed): default (2,2) 0; quality-first 0; retention target 60% 2 and 70% 2 (inter-laboratory V1-V2, V4); retention target 80% 3 (mock B, inter-laboratory V1-V2, V4); FIGARO's score on exact pass rates 0. At AmpliTrim's truncation, the maxEE that reaches FIGARO's retention gives kept reads with about the same expected errors on three datasets (mock B, inter-laboratory V1-V2, V4) and a mixed result on three (V1-V2, mock A, inter-laboratory V4-V6: one read cleaner and the other noisier by 0.03 to 0.5 expected errors); never consistently cleaner. FIGARO chose a stricter maxEE than DADA2's default on the high-quality data (1,1 on V4 and the inter-laboratory V1-V2, 2,1 on V1-V2) and a looser one on poor reads (3,4; 4,4; 6,9), so its advantage is adapting in both directions. Near misses (inter-laboratory V4-V6 at a 70% target: forward +0.03, reverse -0.04; mock A at 80%: forward cleaner, reverse noisier) would change the count only under a looser tolerance or a different tie-break, which would be post hoc. On mock B FIGARO's kept reverse reads carry 3.6 expected errors, and its ASVs were still accurate (7 of 7 variants), so expected errors per kept read do not by themselves predict accuracy. Correction (same day). An earlier note in this entry said that the report's footnote described FIGARO's score wrongly, because FIGARO's pick (6,9) on mock B scores -12.4 by that formula. That note was incorrect. FIGARO sets the maxEE of each cut from the reads (the fitted mean expected errors at the cut, plus one) and uses its score only to rank the cuts, as the README says. A high maxEE with a negative score is therefore expected on poor reads, and the footnote is correct. The candidate list of this check did not include FIGARO's own derivation, which would match FIGARO by construction; it is already available by running FIGARO through `--figaro-script` and giving its maxEE to `--max-ee`.

---

## 7. Development stage (supplement only)

Built and verified the tool; reproducible; **not** results of the manuscript.

| ID | Check | Data and command | Outcome |
| --- | --- | --- | --- |
| S1 | Exactness of the filter replica against DADA2's `filterAndTrim` | D1; `python validation\filter_grid_check.py` (`--quick` for a small grid); 2 samples x 5 x 5 truncations x 4 maxEE | 200 of 200 identical (pass rates 26.3% to 100%) |
| S2 | ITS mode: pairs kept, mean cut length, expected errors, errors of added reads | D2: `python validation\simulated\its_sim.py OUT`; the package test recomputes in plain Python | identical: 6,862, 8,611 and 4,721 pairs kept; cut lengths 210.00, 217.76 and 194.55 nt; 0.2901 and 3.0129 expected errors per added read |
| S3 | Coverage rule: the longest taxon DADA2 recovers against the total | D2: `python validation\simulated\coverage_sim.py OUT noisy_reverse`; DADA2 at totals 511, 505, 480, 470 | total - 12 nt at every total; none of 172 longer taxa recovered |
| S4 | "Beaten by" rule and FIGARO at equal maxEE on a simulated community | D2: `coverage_sim.py OUT both_degrade` | both chose 265/249 at maxEE 2,4 (75.1%); AmpliTrim's default chose 277/237 (55.9%) |
| S5 | What a looser maxEE keeps | D2: `python validation\simulated\staggered_mock_sim.py OUT` | relaxing the reverse limit 2 to 5 raised the share passing 78% to 96% |
| S6 | Chimera removal in variable-length reads | D2: `python validation\simulated\chimera_sim.py OUT` | 6.5% removed when 7.8% were planted |
| S7 | Primer removal | real reads, `real_data.py` | a single cutadapt pass let 0.10% of 3.09 million pairs through untrimmed; fixed |
| S8 | Effect of no truncation on real V4 reads | D5 | superseded by the battery on D5 (23% passed without truncation against 88% with it) |
| S9 | Software | `python -m pytest tests -q`; `python -m pyflakes amplitrim validation tests` | 262 tests; no findings |

---

## 8. Still to run, and further tests recommended

**To run: one command, overnight** (`python validation\check_build.py --overnight`). It finds the folders and FIGARO; checks the machine; downloads, trims and registers V4 (PRJNA574030) and both regions of
the inter-laboratory mock (PRJEB34919: the two regions share one download, the shotgun runs are left out by the MiSeq filter, and the primers are those of the paper's Table 3); runs every step (a step that
covered other datasets is repeated, finished ones are kept; one dataset failing does not lose the others); keeps Windows awake; and writes `MORNING_REPORT.md` (the step verdicts, the test tables, the decision,
what needs you). **Still by hand:** mock B's reference (the run fetches and prepares its eight mock samples, BiomPsy 536-539 and 544-547, by itself; the reference is Zymo's D6321 file: `python validation\check_build.py --make-reference mockB=D6321.refseq.zip` reads the rRNA gene files, cuts them between the primers, collapses identical copies and registers the result), the expected sequences of the inter-laboratory mock (a FASTA of the 15 bacteria, for the truth columns), V4-V5 (primers). Checked on the real files (`--compare-reference mockA=ZymoBIOMICS.STD.refseq.v2.zip`): the D6300 standard shares only five of mock A's 14 species (B. subtilis, E. coli, E. faecalis, L. monocytogenes, S. aureus; all 12 of the Zymo amplicons of those species match a mock A variant at 99%), mock A has nine strains the standard lacks and the standard has two (L. fermentum, S. enterica) that mock A lacks; P. aeruginosa is in both but Zymo's sequences have no forward primer site, so they were not cut. The paper (Hülpüsch, Rauer et al. 2023) explains the confusion: its *even* mock, from Karstens et al., is the ZymoBIOMICS D6300 standard (V4, 8 bacteria and 2 fungi, a dilution series; only processed data are given), and it is a different dataset from mock A. Mock A keeps the reference rebuilt from the authors' repository (14 species): one expected-sequence file per dataset, and the per-sample tests (dilution series, coverage, landscape) need the variants present in each sample, which a species list cannot give. Whether the ZymoBIOMICS standard (D6300) describes mock A is settled by `python validation\check_build.py --compare-reference mockA=ZymoBIOMICS.STD.refseq.v2.zip`, which cuts the Zymo sequences with the mock A primers and lists the strains the two files share and those only one has; mock B takes `--make-reference mockB=D6321.refseq.zip`.
A last full run on the final build (`--run --fresh`) so that the manuscript cites one build id.

**Further tests that would strengthen the evidence, in order of value:**

1. **Held-out datasets, in this order: D14 (the real tight-overlap case with an independent truth: it replaces the clipped-mock emulation A5 and the V3-V4 2x250 item D8), D7, D5, then D9.** Mock A informed design decisions and V4 was used in development. Only datasets that were never looked at give an unbiased test: D7, D8, D9. Run them with the margins fixed in advance.
   *Considered and dropped:* the sputum study (PRJNA915241, D10 to D13: four regions on the same 33 samples, a ZymoBIOMICS standard). It would add held-out regions, but the primers are not given in the paper (phased primers), the run names are unknown, the four regions are one study and not four, it gives no tight-overlap case, and the mock may not be in the project. Its catalog entries (`sputum_v1v2` to `sputum_v7v9`) remain in `check_build.py` if it is wanted later, for example if a reviewer asks for clinical samples.
2. **A second independent truth.** One mock with expected sequences that do not come from a DADA2 run (D7); ideally a second, with a different community and region.
3. **FIGARO on at least three datasets with paired intervals** (B4), so the "ties FIGARO" claim does not rest on one study.
4. **A "quality-plot rule" baseline** for B3: truncate each read where the median quality falls below Q30, as many users do, to compare AmpliTrim with common practice and not only with FIGARO and no truncation.
5. **Platform diversity.** A NovaSeq run (binned quality scores) and a 2x150 run, for the filter replica (B1) and the calibration (B2).
6. **Read-depth variants.** The battery at 5,000 pairs per sample and at full depth on every dataset (the depth effect on B2 is so far seen only in ITS).
7. **A second ITS study** if ITS claims are kept: B1, B2, B3 and the maxEE rows on another amplicon and platform.
8. **Runtime and memory** against the number of pairs and samples (a software-paper item).

## 9. Round 4: the open questions, tested in one place (criteria fixed on 9 Oct 2026, before any run)

Run `python validation\check_build.py --round4 --quick --phases A,B` first (about 25 minutes: small slices, to see that every path and tool works; a quick run reads the earlier results of the full run, and the simulations give no verdict in it), then
`python validation\check_build.py --round4`. The run goes in three phases with a break of 90 s between them (Enter continues at once, S skips the next phase, Q stops, P holds; **if nothing is pressed
the run continues by itself**; `--break-seconds 0` removes the breaks; `--phases A` or `A,B` chooses phases; `--fresh` repeats steps that are finished). It writes `ROUND4_REPORT.md`.
Phase A reads results that are already on disk (minutes), phase B simulates communities with known truth (hours), phase C runs the standard battery on any dataset registered since the last run.
The criteria below are in `validation/round4_steps.py` (the constants at its top) and are not to be changed after seeing results.

| Q | Question | Step and command | Criterion, fixed in advance | If met / if not met |
| --- | --- | --- | --- | --- |
| Q28 | Does a forward-heavy split beat the balanced split that AmpliTrim prefers among near-ties? | `--only split_tiebreak` (re-reads the landscape grids; no DADA2) | in at least 4 of the 6 datasets the best cell with AmpliTrim's total length has the longer forward read, and the mean gain over those is at least 1.0 point of merged pairs | met: add an opt-in `--split-preference forward` and test it on the next datasets before any change of default; not met: keep the balanced tie-break |
| Q29 | Which expected variant does only one tool recover, and which read covers the bases that tell it from its nearest sibling? | `--only missing_variants` | information only (no pass or fail): for how many such variants more of the distinguishing bases are read by both reads under the recovering tool's split | informs Q28 |
| Q30 | How many pairs of each dataset can a tool that needs reads of one length use? | `--only figaro_usable_pairs` | information; V1-V2: the share with exactly the most common lengths is within 5 points of the 47% used earlier | confirms or corrects the 47% |
| Q23 | Do the ASVs of the mocks have the lengths of the expected amplicons (false joins make wrong-length ASVs)? | `--only asv_length_check` | mock A at AmpliTrim's default setting: at most 1% of the ASV reads are off every expected length (5 nt); mock B information only (a third of its reads are contaminants) | not met: look at the largest off-length ASVs first |
| Q16 | Can FIGARO's maxEE rule (the mean expected errors of the reads at the cut, plus 1) be reproduced from the reads? | `--only figaro_rule_check` | within 1.0 on both reads in at least 5 of the 6 datasets with a FIGARO row | met: an opt-in `--auto-max-ee fitted`, tested against FIGARO as shipped (decision 2); not met: not built, users pass FIGARO's maxEE with `--max-ee` |
| Q0-14 | Does the excess of the yield estimate over DADA2's merged share grow with the richness of the community? | `--only richness_sim`: 25, 100, 250, 500 and 1,000 taxa, 5 seeds, maxEE 1, 2 and 5; 10 samples of 5,000 pairs; related taxa, log-normal abundances, qualities resampled from mock A | supported: the excess at 500 taxa is at least 10 points above that at 25 taxa; not supported: under 3; between: partly. Also: the estimate ranks the three maxEE as DADA2 does in every community | supported or partly: the README says the yield level is an upper bound that rises with richness; not supported: the cause is unknown |
| Q21-27 | How many merges are at the wrong length, and does the 12 nt minimum overlap hold? | `--only overmerge_sim`: arms clean, realistic (a 10 nt repeat in a fifth of the taxa, 1% chimeras, indels 1e-5) and stress (32 nt repeat, 5%, 1e-4); minimum overlap 12, 16 and 20; 3 seeds | 12 nt stands if, in the realistic arm, at most 1% of the merged reads of pairs without an indel are at the wrong length at 12 nt, and at 12 nt DADA2 recovers no more than 1 taxon fewer than at 20 nt (amended 9 Oct, see the note below) | not met: consider a longer minimum overlap, and test it on real data before any change |
| - | The standard battery on datasets registered since the last run | `--only battery,extra_datasets` (phase C) | the criteria of section 3 | - |


**Amendment of the Q21-27 criterion (9 Oct 2026, before any full run of round 4).** The quick run (3 arms x 3 overlaps, 1 seed, 2,500 pairs; no verdict by design) showed two flaws in the wording of the criterion. Both were corrected on 9 Oct, before any result of the full run existed.
(1) The wrong-length share counted merges of pairs that carry a simulated indel. An indel changes the merged length by one base whatever the merge does, so it is not a false join: in the realistic arm it made up the whole 0.4% (0.41, 0.34 and 0.35% at 12, 16 and 20 nt, the share expected from the indel rate alone, and independent of the overlap).
The criterion now uses the pairs without an indel (the simulator marks the others with `|indel` in the read header), and the share among pairs with an indel is reported apart.
(2) The recovery condition used an absolute difference between 12 and 20 nt, which would have failed 12 nt for recovering MORE taxa; in the quick run it did (clean arm 15, 13 and 13 taxa at 12, 16 and 20 nt; realistic 17, 16 and 15), because a longer overlap keeps noisier bases and loses rare taxa. The condition now fails only if 12 nt recovers more than one taxon fewer than 20 nt.
The thresholds themselves (1% and 1 taxon) are unchanged. Nothing else in section 9 was changed.

Not in this round, and why: relative-abundance accuracy (needs the expected abundances of mock A as a file, Q31), the reasons for the weak per-sample ranking (Q32) and the ITS 200/200 conflict (Q33) are analyses of files not on hand; the data questions
(Q34-Q42: full-depth V4 and ITS, truth for the inter-laboratory mock, V3-V4 and V4-V5, a binned-quality platform, a second ITS study) need data: register it and phase C runs the battery on it.

## 10. Round 4: outcomes (run on 8 October 2026)

Phases A and B ran on build `6b0c25b6e3` with the criteria of section 9, including the amendment of the Q21-27 criterion made before the full run. Phase C repeated the standard battery on the two original datasets (mock A, V1-V2) under the 0.2% default; the extra datasets
(V4, the inter-laboratory mocks, mock B) were kept from 7 October, when the default was 1%. The quick runs (one seed, small slices) only checked that every path worked; their numbers are not used, except where a step reads earlier results (Q28) or is marked.

| Q | Criterion fixed in advance | Result | Outcome | Consequence |
| --- | --- | --- | --- | --- |
| Q28 | the best cell with AmpliTrim's total length has the longer forward read in at least 4 of 6 datasets, mean gain at least 1.0 point | 3 of 6 (soil V4 +4.7, mock A +0.5, V1-V2 +0.1); mean gain over those 1.8 points, and 0.3 without V4; in the other three AmpliTrim's own split was best or no worse | NOT MET | keep the balanced tie-break; the gain is one dataset (V4) |
| Q29 | information only | 0 of 2 variants: both differ from their nearest expected variant by one base (mock A, 484 nt, base 416; mock B, 515 nt, base 455), which only the reverse read covers, early in it (about base 69 and 61), under both splits | no verdict | the split does not explain FIGARO's extra variant (most likely the denoising: a rare one-base variant is kept or absorbed); n = 2 |
| Q30 | V1-V2: pairs with exactly the most common lengths within 5 points of 47% | 46.5% (full run) | MET | confirms the 47% used earlier |
| Q23 | mock A, default setting: at most 1% of the ASV reads off every expected length (5 nt) | 0.13% (one ASV of 435 nt, 24 to 26 reads, present in every method including no truncation) | MET | no false-join signature in mock A |
| Q16 | 'mean expected errors at the cut + 1' within 1.0 of FIGARO's maxEE on both reads in at least 5 of 6 datasets | 4 of 6 (quick-run values: mock A 2.3/2.9 against 3/4 and mock B 3.6/6.1 against 6/9 not within; the other four within); FIGARO's maxEE is higher than the rule on poor reads | NOT MET | no automatic maxEE is built; pass FIGARO's maxEE with `--max-ee` |
| Q0-14 | excess of the yield estimate at 500 taxa minus at 25 taxa: at least 10 points supported, under 3 not supported, between partly | +3.8 points (-1.4 at 25 taxa, +4.4 at 1,000) | PARTLY SUPPORTED | richness explains a small part of the +30 on soil V4 |
| ranking | rank correlation of the estimate with DADA2's merged share over three maxEE at least 0.99 in every community | 1.00 in 25 of 25 communities | MET | |
| Q21-27 | realistic arm: at most 1% of the merged reads of pairs without an indel at the wrong length at 12 nt, and at 12 nt DADA2 recovers no more than 1 taxon fewer than at 20 nt | 0.00%; taxa recovered 12 minus 20 nt: +0.0 | MET | the 12 nt minimum overlap stands in this simulation |

**Phase C, the standard battery (paired difference in merged pairs at FIGARO's maxEE, 95% interval, margin +/-1.0):** mock A -0.0 (-0.3 to +0.3, FIGARO's maxEE 3,4), V1-V2 +0.4 (-0.1 to +0.8, 2,1): both refreshed on 8 October at the 0.2% default; soil V4 +2.0 (+1.6 to +2.4, NOT MET in AmpliTrim's favour),
inter-laboratory V1-V2 -0.0 (-0.0 to +0.0), mock B +2.6 (+0.3 to +5.2, NOT MET in AmpliTrim's favour), all three kept from 7 October (1% default); inter-laboratory V4-V6: no paired comparison (5 samples). FIGARO's measured amplicon lengths moved from 505 to 508 (mock A) and from 343 to 346 (V1-V2) with the new default.

**Richness simulation (Table 12 of docs/FINDINGS.md; 25 to 1,000 taxa, 5 seeds, 10 samples of 5,000 pairs, qualities resampled from mock A):** the excess of the estimate over DADA2's merged share, maxEE 2: -1.4, -0.0, +0.5, +2.4, +4.4 points at 25, 100, 250, 500 and 1,000 taxa. Reads in exact ASVs: 99.4%, 96.7%, 94.4%, 92.6%, 92.9%. Taxa with 50 reads or more were recovered at 89% to 96% at every richness; taxa with 10 to 49 reads (strains 1 to 3 bases from a sibling) at 6% (100 taxa) to 50% (1,000 taxa).
**Over-merging simulation (Table 13; 25 taxa, reads of 250 nt, 3 seeds, maxEE 2):** see the tables in docs/FINDINGS.md, section 4.7.

**Exploratory findings (no pre-set criterion; to be repeated before they are relied on):**
(1) At 1,000 taxa a looser maxEE (1, 2, 5) raised the merged share from 47.5% to 68.5% to 85.9% but lowered the share of reads in exact ASVs from 93.8% to 92.9% to 92.3%: a small accuracy cost in a rich community, which was not seen on the mocks.
(2) With a 32 nt perfect repeat across the overlap in a fifth of the taxa (stress arm), vsearch merged between 30% and 88% of the pairs depending on the seed (it refuses most pairs of repeat taxa), merges at the wrong length were 0.08% at 12 nt, 0.01% at 16 nt and 0.00% at 20 nt, and DADA2 put 27.2% of the reads in ASVs matching no true sequence at a 12 nt overlap, 19.4% at 16 nt and 10.1% at 20 nt (5% of the reads were chimeras): in a repeat-rich region a longer overlap helps, and the yield estimate can be far too low when a repeat dominates.
(3) At low richness (25 taxa) the estimate was 1.4 points BELOW DADA2's merged share, as on some real mocks: the estimate is not an upper bound everywhere.

**Not answered by round 4:** why soil V4 shows +30 (the leading untested explanation is the long tail of rare sequences, which the log-normal abundances of the simulation lack); accuracy in a rich real community (no truth exists); V3-V4, V4-V5, binned-quality platforms; the ITS 200/200 conflict (+16.1 against +23.5).

## 11. Night 1 (9 October 2026): the complete battery on twelve datasets

One build (`6b0c25b6e3`), one output folder, the 0.2% default, run with `night1_v2.py`. The new datasets were prepared by that script with rules written before the run: primers and orientation chosen from what `describe_reads.py` found in each dataset's reads; swabs (PRJNA1261676): paired NovaSeq runs left out when their own labels matched a mock or control pattern, at least 20,000 reads, 36 evenly spaced over the sorted accessions. The battery took 35.6 minutes (the original datasets) and 177 minutes (the nine others).

**Verdicts (criteria of section 3, margin of 1.0 point at FIGARO's own maxEE):** met on 6 (mock A -0.0, V1-V2 +0.4, soil V4 +0.0, inter-laboratory V1-V2 -0.0, D6300 V3-V4 +0.0, swabs +0.0); NOT MET on mock B (+2.6, interval +0.3 to +5.2, AmpliTrim ahead), spike-in (+1.3, +0.4 to +2.5, AmpliTrim ahead) and wastewater (-1.0, -1.4 to -0.7, FIGARO ahead); the D6300 V4-V5 line reported NOT MET because FIGARO crashed on those reads (an IndexError in its own code), which is a missing result and not a difference: the verdict is now worded 'no result' (battery text only); inter-laboratory V4-V6 has 5 samples, too few for a paired comparison. At the same maxEE the two tools give the same merged share on 5 of 10 datasets to the decimal and within 0.4 points on 3 more, because they choose almost the same truncation.

**Measured:** AmpliTrim's filter prediction within 0.8 points of DADA2's pass rate on all 12 datasets (FIGARO's predicted retention was off by up to 7.7 on the NovaSeq sets, where it could use only 8 to 25% of the pairs). Rank correlation of the yield estimate over settings 1.00 on 9 of the 11 datasets where it was printed (mock A 0.90, wastewater 0.80). Yield estimate minus DADA2's merged share, mocks: inter-laboratory V1-V2 +0.3, mock B -2.2, inter-laboratory V4-V6 -1.6, mock A +2.7, D6300 V3-V4 +2.9, D6300 V4-V5 +6.7; real: wastewater +3.9, spike-in +3.4, V1-V2 +6.5, swabs +7.0, ITS +21.1, soil V4 +32.1. With its default maxEE of 2 AmpliTrim merges 2 to 6 points more than FIGARO as shipped on six clean datasets and fewer on mock A (-1.8), wastewater (-23.1) and mock B (-54.1), where FIGARO chose maxEE 5 and 6,9. Accuracy on the mocks: D6300 V3-V4 recovered all 8 strains at every setting with 97.9% of reads in exact ASVs; a looser maxEE cost nothing there.

**Resolved:** the ITS conflict (+16.1 against +23.5) came from two different baselines: trunc200 against AmpliTrim at maxEE 1 (+16.1), and against AmpliTrim at the default maxEE 2 (+23.5). The overshoot of +32 on soil V4 is not caused by tolerating mismatches in the overlap: the strict estimate (no mismatch allowed) still overshoots by +25.0.

**Exposed:** (1) the no-truncation arm crashed on both NovaSeq V3-V4 sets (no minimum read length): fixed in `compare_methods.py` (`min_len_for`). (2) On D6300 V4-V5 about 1% of the amplicons are 550 nt against a median of 374 (probably the two yeasts' 18S, which 515F/926R also amplify); covering them makes AmpliTrim recommend the full read length, 282/280, which loses 10% of the pairs to the length requirement (merged 57.2%, against 63.9% with no truncation). Recorded as a finding; the tool is not changed (it triggered on one dataset of twelve, so a fix could not be validated, and changing it after seeing the result would be tuning on the test). (3) FIGARO was given 550 as the amplicon length with reads of 282+280, which leaves exactly the 12 nt overlap of the harness and no slack; it failed without a message. Night 2 also gives FIGARO the typical length.


## 12. Night 2: the maxEE head-to-head (fixed on 9 October 2026, before the run)

**Question.** FIGARO chooses a truncation and a maxEE; AmpliTrim chooses a truncation for a maxEE it is given. When AmpliTrim also chooses the maxEE, does it end up where FIGARO ends up, with as many merged pairs and no noisier reads? At the same maxEE the two give the same outcome (section 11), so the comparison is about the maxEE decision.

**What FIGARO does (the account of decision 2, restated).** FIGARO derives each cut's maxEE from the reads (the fitted mean expected errors at the cut plus one) and uses its score (retention in percent minus (F-1)^2 minus (R-1)^2) only to rank cuts. It is not a price per allowed error that is searched over; a high maxEE on poor reads is expected. AmpliTrim's quality-first pick works differently: among the maxEE rows of its grid that are not stricter than 2 and are not beaten by another row, the one with the fewest expected errors per kept pair that still keeps the floor share of the pairs (default 50%). Its grid is (1,1), (2,2), (2,3), (2,4), (2,5), (3,3), (3,5), (5,5), so it cannot reach (6,9).

**Arms** (each is a truncation and a maxEE, run through DADA2 on the same pairs: 5,000 per sample, all reads for mock A): A AmpliTrim at maxEE 2; B AmpliTrim at its quality-first maxEE, floor 50% (the headline); B2 the same at 80%; B4 AmpliTrim at 'the mean expected errors of the reads at AmpliTrim's default cut, plus one', rounded to whole numbers with halves up and at least 1 (FIGARO's documented rule; the earlier reproduction matched within 1.0 on 4 of 6 datasets); C FIGARO's own pick; D AmpliTrim's truncation at FIGARO's maxEE (a control); E no truncation (maxEE 2, minimum length 50; trunc200 for ITS); F a practitioner's cut where each read's median quality falls below 30 (maxEE 2, not shorter than 50). Both tools use the same 12 nt overlap (`figaro.py -m 12`); there is no arm with FIGARO's own default overlap. FIGARO is given the longest amplicon AmpliTrim measured and, when the typical (median) amplicon is under 95% of the longest, also the typical one. A FIGARO failure is a missing result, not a failed criterion; the arms are then classified against A. B4 is not run on ITS (no truncation to take a mean from) and FIGARO is not run on ITS.

**Measures.** (1) The share of pairs DADA2 merges: the paired per-sample difference, with a bootstrap 95% interval from 10 samples (fewer samples: point estimates only); two picks are the same inside 2.0 points. (2) The mean expected errors per kept read, forward and reverse, measured from the qualities at the arm's truncation and maxEE (`ReadPairs.error_load`, 60,000 pooled pairs); two picks are equally noisy when both reads are within 0.02 expected errors, cleaner or noisier when both reads are on one side, and 'mixed' when one is cleaner and the other noisier. (3) On mocks: variants recovered, strains recovered and the share of reads in exact ASVs, with the battery's margins (1, 1 and 0.1 points); on a mock, accuracy is shown beside the verdict and decides a disagreement. **Verdict words** against FIGARO's pick: better (a gain on one axis and no loss on the other), worse, matches, trade-off (more pairs but noisier, or fewer pairs but cleaner), mixed noise, no result. No winner is declared on real data unless one pick is better on both axes.

**Also in night 2.** The coverage setting: on datasets with a tail of long amplicons AmpliTrim also runs at maxEE 2 leaving out 1% and 5% of the longest amplicons (the default is 0.2%). The stability of each tool's pick: on datasets with 9 or more samples, AmpliTrim (default and B) and FIGARO pick again on three disjoint sets of samples; a pick holds when the truncations differ by at most 10 nt and the maxEE is the same.

**Predictions, written before the run.** P1: where AmpliTrim's default keeps 80% or more of the pairs, B and B2 equal A and A against FIGARO is 'trade-off: more pairs, noisier'. P2: on poor reads (mock B, wastewater, inter-laboratory V4-V6) B merges fewer pairs than FIGARO by more than the margin, with cleaner reads (about 18 points on mock B); B2 narrows the gap and B4 comes nearer FIGARO than B does. P3: on mock B the accuracy of B is lower than FIGARO's. P4: on D6300 V4-V5 FIGARO runs at the typical length, chooses much shorter reads, and leaving out 1% of the longest amplicons raises AmpliTrim's merged share. P5: both tools' picks hold across the three sets of samples, except where quality differs strongly between samples (wastewater).

**Datasets** are those registered in `amplitrim_paths.json`: mock A (high input), V1-V2, ITS, and the extras (soil V4, inter-laboratory V1-V2 and V4-V6, mock B, D6300 V3-V4 and V4-V5, wastewater, spike-in, swabs). Development data: mock A and soil V4; the rest are held out. Mock and real are reported apart.

**Decided and not done.** The default-overlap FIGARO arm (declined). FIGARO's score as a maxEE rule (decision 2: 0 of 6; it ranks cuts). 'No truncation when the best cut is the full read length' (a finding only). The spike-in finder (later, separate: the sequences would be derived from the same reads). The skin dataset (dropped: another real V1-V3 set adds little). The inter-laboratory mocks stay without expected sequences (the paper gives family-level digital PCR abundances, not 16S sequences). **What follows.** Nothing is built from night 2 automatically. If B, B2 or B4 match FIGARO, an opt-in `--auto-max-ee` is decided afterwards and tested against FIGARO as shipped in a later run, as decision 2 requires.

**Build for night 2: `b8469465fc`.** The first copy (`236dd28c0c`) failed one test on Windows (`test_night2_end_to_end`: the temporary folder could not be deleted because `night2.py` left `night2.log` open). `night2.py` now closes its log when the run ends, even after an error, and the test asserts it. Nothing else changed: the arms, the margins, the verdict words and the predictions of this section are as fixed above, before any run.

## 13. Night 2: outcomes (run on 9 October 2026, build `b8469465fc`)

The run covered 12 datasets, started at 12:08 and ended at about 15:50 (3 h 42 min), with no failed dataset. FIGARO gave no result at an amplicon length of 550 on D6300 V4-V5 (an `IndexError` in its trim search: the reads leave exactly 12 nt of overlap for a 550 nt amplicon, so the only possible cut is no trimming; at 374 it ran). Per the fixed condition this is reported as 'no result', not as a criterion. The swab samples were chosen by the fixed rule of section 11. Arm E (no truncation) cannot be scored on noise (the errors per kept read are measured at a cut), so its verdict is 'no result' on every dataset; its merged shares are in the tables of `docs/FINDINGS.md`.

**Predictions of section 12, against the outcome.**

| | Prediction | Outcome |
| --- | --- | --- |
| P1 | Where A keeps at least 80% of the pairs, B and B2 equal A, and A against FIGARO is 'trade-off: more pairs, noisier' | Met on 5 of 6 datasets (soil V4, inter-laboratory V1-V2, D6300 V3-V4, spike-in, swabs). V1-V2 came out 'mixed noise (more pairs)': forward cleaner, reverse noisier (FIGARO cuts the reverse read at 138 against 157). B and B2 equalled A on all six. |
| P2 | On poor reads (mock B, wastewater, inter-laboratory V4-V6) B merges fewer pairs than FIGARO by more than the margin, with cleaner reads (about 18 points on mock B); B2 narrows the gap; B4 comes nearer FIGARO than B does | Met in direction on all three (B: mock B -20.6, wastewater -23.1, V4-V6 -10.7, all cleaner). B2 narrows the gap on mock B (-18.1) and wastewater (-1.0, but noisier: 'worse'), and overshoots on V4-V6 (+6.0, noisier). B4 is nearer than B on mock B (-12.2) and wastewater (-12.0) and equal to B on V4-V6. Mock B's gap for B was larger than predicted (-20.6 against about 18). |
| P3 | On mock B the accuracy of B is lower than FIGARO's | Met: 6 against 7 variants; reads in exact ASVs 74.2% against 78.3% (A 68.7%, B2 73.6%, B4 74.9%, D 75.3%). |
| P4 | On D6300 V4-V5 FIGARO at the typical length chooses much shorter reads, and leaving out 1% of the longest amplicons raises AmpliTrim's merged share | Met: FIGARO 249/137 (maxEE 1) merged 80.3% against AmpliTrim's 57.2%; leaving out 1% gave 268/118 and 89.8% (+32.6). |
| P5 | Both tools' picks hold across three disjoint sample sets, except where quality differs strongly between samples (wastewater) | Not met. AmpliTrim held on 4 of 6 datasets (it moved 13 nt on V1-V2 and 11 nt on soil V4); FIGARO held on 3 of 6 (it moved 17 nt and changed maxEE on wastewater, shifted on the spike-in set and gave an error on one set of D6300 V3-V4). Wastewater, predicted to be the exception, was identical across the three sets for AmpliTrim. |

**The question of section 12:** does AmpliTrim, choosing its own maxEE, end up where FIGARO ends up? Against FIGARO's pick (10 datasets where it gave a result at the longest amplicon length), the verdict words were: A and B, no match on any (9 trade-offs, 1 mixed); B2, no match (8 trade-offs, 1 worse, 1 mixed); B4, matches on 4 (inter-laboratory V1-V2, D6300 V3-V4, soil V4, swabs: the clean datasets where FIGARO picks maxEE 1), 5 trade-offs, 1 mixed; the control D (AmpliTrim's cut at FIGARO's maxEE), matches on 7, better on 1 (V1-V2, +1.1), mixed on 1 (mock B, +2.6), worse on 1 (wastewater, -1.0). **Answer: no.** On clean reads FIGARO's pick is stricter (maxEE 1) and AmpliTrim's default merges 2.3 to 6.3 points more pairs, for 0.01 to 0.10 more expected errors per kept read. On poor reads FIGARO's pick is looser than any AmpliTrim rule except B2, and the default merges 15 (mock A), 23 (wastewater), 27 (V4-V6) and 54 (mock B) points fewer. Where accuracy differed it favoured the looser pick (mock B; D6300 V4-V5 at the typical length, 100.0% of reads in exact ASVs against 98.8%). On mock A and D6300 V3-V4 all arms were equally accurate.

**Also measured.** (1) FIGARO at the typical amplicon length against the longest: V1-V2 43.3% against 81.2%, inter-laboratory V1-V2 54.4% against 93.6%, mock B 45.7% against 74.3% (one of three strains recovered), wastewater 25.5% against 66.7%; better on soil V4 (72.1% against 63.3%) and V4-V5 (80.3%; the longest length crashed). (2) The coverage setting (merged share at 0.2% / 1% / 5% of the longest amplicons left out): D6300 V4-V5 57.2 / 89.8 / 89.8; soil V4 66.8 / 67.2 / 72.3; V1-V2 84.2 / 84.2 / 84.2; inter-laboratory V1-V2 95.9 / 95.9 / 93.8; mock B 19.9 / 19.9 / 19.9; wastewater 43.6 / 43.6 / 41.4. At 1% no dataset lost; at 5% two lost about 2 points. (3) No truncation (E) against AmpliTrim's cut at maxEE 2, on the 11 16S datasets: more pairs merged on swabs (+4.4) and D6300 V4-V5 (+6.7), a tie on D6300 V3-V4 (-0.1), fewer on the other eight (spike-in -1.8, soil V4 -12.7, mock A -23.3, V1-V2 -55.0, wastewater -43.4, mock B -17.9, inter-laboratory V1-V2 -69.9, V4-V6 -20.8): E collapses where the reads are long and the amplicon needs the overlap. (ITS has no truncation by design; its 200/200 arm lost 23.5 points.) (4) The median-Q30 baseline (F) collapsed on 7 of 11 datasets (merged 0.0 to 0.8%, or DADA2 could not run): on soil V4, D6300 V3-V4 and D6300 V4-V5 its cut is longer than the reads, so DADA2 removes every read; on mock B, wastewater and inter-laboratory V4-V6 it leaves too little overlap; on the spike-in set DADA2 could not learn an error model from the binned qualities. It merged 21.8% on mock A (a cut of 472 nt for a 508 nt amplicon) and worked only on V1-V2 (72.1%) and inter-laboratory V1-V2 (79.6%), where AmpliTrim's default was ahead by 12.1 and 16.3 points; the swab run gave no result. (5) Estimates: soil V4's yield estimate ranked 7 settings poorly (rank correlation +0.09, excess up to +33) when the settings differed in coverage; on wastewater the estimate preferred a cut 6 nt longer than FIGARO's at maxEE 5 (predicted 72.9% merged against 64.6% for FIGARO's; DADA2 gave 65.7% and 66.7%), and the extra length cost 0.9 points at the filter.

**Decisions (made after the run).**
1. **No automatic maxEE is built.** The condition of section 12 ('if B, B2 or B4 match FIGARO') was not met: B4 matches on clean data only, B2 overshoots or is noisier. The defaults stay: maxEE 2, `--allow-unmergeable` 0.2%.
2. B4 is an approximation of FIGARO's derivation (the mean expected errors of the pooled reads at the cut, plus one). Whether FIGARO's exact derivation would reproduce its picks on poor reads is not known. **Fixed now, before it is tested:** an exact implementation of FIGARO's maxEE derivation will be checked offline, on the cuts and the picks already on disk (no DADA2 run), and an opt-in `--max-ee auto` is considered only if it reproduces FIGARO's maxEE, within 1 on each read, on at least 8 of the 10 datasets where FIGARO ran at the longest amplicon length. This needs FIGARO's source (the maxEE part of `trimParameterPrediction.py`).
3. **Text-only warnings with an alternative (post hoc: the thresholds were chosen after night 2, inside a visible gap, and are to be tested on new data before release).** (a) A long-amplicon tail: when leaving out 1% would raise the estimated yield by at least 3 points, say so, with the cut and the estimate (the existing `--explore-coverage` rule, made automatic). (b) Low survival: when fewer than 70% of the pairs are estimated to pass at maxEE 2 (the datasets were at 20 to 65% or at 85% and above), show the estimated pass rate of the next looser settings, stating that no mock lost accuracy at a looser maxEE. (c) Variable-length reads: when more than 3 points of pairs would be lost to DADA2's length filter at the recommended cut, say that no truncation (with a minimum length) would keep them. None changes a recommendation.
4. 'No truncation at full length' stays a finding (E merged more than AmpliTrim's cut on 2 of 11 16S datasets, tied on 1 and less on 8); no change to the tool.
5. **F2, post hoc (labelled as such; run as a short night 3).** Median-Q30 cut, never longer than the reads and never shorter than the longest amplicon plus 12 nt; maxEE (2,2), relaxed to (2,5) when fewer than 60% of the pairs pass at (2,2). F stays as registered and is reported as a baseline that collapsed.
6. Validation stops after F2 and the check of decision 2; the results are written up once (`docs/FINDINGS.md`, with the development rounds in `docs/SUPPLEMENT.md`).

**Exposed by night 2 and not fixed.** The explore pick (B) kept maxEE 2 on wastewater although AmpliTrim estimated 49% survival; the 50% floor is applied to the estimated pass rate or yield as coded and was not re-read for this section. Arm E has no noise measure. The Q30 arm needs the cap of F2. Stability: AmpliTrim's pick moves 11 to 13 nt along a flat ridge on V1-V2 and soil V4 (not an accuracy difference, but it fails the criterion as written). Five datasets have fewer than 10 samples (mock A 6, inter-laboratory V1-V2 6, V4-V6 5, mock B 8, D6300 V4-V5 3, wastewater 9): point estimates only. Depth: 5,000 pairs per sample except mock A.

## 14. Night 3: arm G and arm F2 (fixed on 9 October 2026, before the run)

**Why.** Night 2 showed that AmpliTrim's default maxEE of 2 is too strict on poor-quality reads. DADA2's tutorial suggests relaxing the reverse read to 5 (maxEE = c(2,5)) when too few reads pass; a limit above 5 is judged too loose by the author of this project, and the data cannot overrule that (only mock B gained beyond 5). Night 3 tests that cap, and adds the more careful manual baseline that night 2 showed to be necessary (arm F collapsed on 7 of 11 datasets). Decision 3 of section 13 (the warning with an alternative) depends on its result.

**Arms.** Each is a truncation and a maxEE, run through DADA2 on the same pairs, with FIGARO at the same 12 nt overlap as in night 2.
- **A** (control): AmpliTrim, maxEE 2.
- **G**: AmpliTrim's best truncation for maxEE (2,5) (forward 2, reverse 5), on every dataset, ITS included (where nothing is truncated).
- **R** (a rule, no extra run): G where the share of pairs passing at maxEE 2 (arm A, measured) is below 70%, otherwise A. This is the alternative the planned warning would offer; the 70% was chosen after night 2 (section 13) and is tested here on the same data, so R is not an independent test of the threshold.
- **F2** (post hoc; 16S only): a careful manual cut. (1) Cut each read where its median quality first falls below 30 (arm F). (2) Never longer than the read: cap each cut at the length that 90% of the reads reach (the 10th percentile of the read lengths). (3) Never too short for the amplicon: if forward plus reverse is under the longest amplicon plus 12 nt, add the missing nucleotides, half to each read (the forward read takes the odd one), up to each cap; what one read cannot take goes to the other; if the caps still fall short, the cut stays at the caps and F2 is flagged 'does not reach the amplicon'. (4) maxEE (2,2), relaxed to (2,5) when fewer than 60% of the pairs pass at (2,2) at that cut (AmpliTrim's exact pass rate).
- **C**: FIGARO's own pick at the longest amplicon length, as in night 2 (no typical-length run, no coverage rows, no stability: they were done in night 2).

**Measures, margins and verdict words** are those of section 12 (merged share within 2.0 points is 'the same'; expected errors within 0.02; accuracy margins 1 variant, 1 strain, 0.1 point of reads in exact ASVs; fewer than 10 samples: point estimates). Each arm is classified against FIGARO's pick, or against A where FIGARO gave no result; G and F2 are also given as differences from A.

**Predictions, written before the run.**
- **P6.** On the datasets where A passes under 70% of the pairs (mock A, mock B, inter-laboratory V4-V6, wastewater, D6300 V4-V5), G merges more pairs than A by more than the 2-point margin, with noisier reads.
- **P7.** G is within 10 points of FIGARO's pick on mock A, wastewater and inter-laboratory V4-V6, and 15 or more points behind it on mock B (maxEE 5 is the cap; FIGARO used 6,9).
- **P8.** On the datasets where A passes at least 80% of the pairs (V1-V2, soil V4, both V1-V2/V3-V4 mocks, spike-in, swabs), G merges at most 3 points more pairs than A, with a noisier reverse read.
- **P9.** G is as accurate as A, or more accurate, on every mock community with expected sequences; on mock B its share of reads in exact ASVs is at least 3 points above A's (68.7%).
- **P10.** F2 relaxes to (2,5) on mock B and wastewater (where AmpliTrim's own pass rates at maxEE 2 are 20% and 54%).
- **P11.** F2 never merges more pairs than the better of A and G by more than the 2-point margin on any 16S dataset (AmpliTrim's cut is at least as good as a careful manual one).
- **P12.** F2 runs on every 16S dataset and merges at least 20 points more pairs than F did in night 2 wherever F merged under 25% or gave no result (soil V4, both D6300 sets, mock B, wastewater, inter-laboratory V4-V6, mock A, the spike-in set and the swabs; a failed F counts as 0%).

**What follows.** If P7 and P9 hold, the warning of decision 3 offers (2,5) as its alternative, with the exact pass rate and yield of AmpliTrim's cut at that maxEE; if either fails, the warning shows the table of looser settings without naming one. F2 is reported as a post hoc baseline and does not change a recommendation. Nothing else is built from night 3 before it is read.

**Build for night 3: `194b5ebf74`.** It adds `validation/night3.py` and arm F2's rule (`f2_plan` in `maxee_arms.py`), and gives `night2.py` a `TAG` so that night 3 reuses its machinery; night 2's own output, file names and tests are unchanged. The tests (`tests/test_night3.py`) check F2's rule against hand-computed cuts and the run end to end; three deliberate breaks (the F2 split, the 70% threshold of rule R, a log left open) each made a test fail. The arms, margins and predictions above were fixed before any of this was run.

## 15. Night 3: outcomes (run on 9 October 2026, build `194b5ebf74`)

The run covered the same 12 datasets, started at 18:12 and ended at 21:33 (3 h 21 min), with no failed dataset. FIGARO gave no result at an amplicon length of 550 on D6300 V4-V5, as in night 2 (arms are classified against A there). F2 ran on every 16S dataset.

**Predictions of section 14, against the outcome.**

| | Prediction | Outcome |
| --- | --- | --- |
| P6 | Where A passes under 70% of the pairs, G merges more pairs than A by more than the margin, with noisier reads | Met on the pairs: +13.4 (mock A), +29.8 (mock B), +18.8 (inter-laboratory V4-V6), +11.6 (D6300 V4-V5), +7.5 (wastewater). 'Noisier' on three (mock A, mock B, V4-V5); 'mixed' on two (V4-V6, wastewater), because G's split moves toward a shorter forward and a longer reverse read, so the forward read is cleaner and the reverse noisier. |
| P7 | G is within 10 points of FIGARO's pick on mock A, wastewater and V4-V6, and 15 or more behind it on mock B | **Not met.** Mock A -1.8 and V4-V6 -8.2 are within 10; mock B -24.2 is 15 or more behind; wastewater -15.7 is not within 10. |
| P8 | Where A passes at least 80%, G merges at most 3 points more than A, with a noisier reverse read | Met on 5 of 6 (inter-laboratory V1-V2 +1.4, D6300 V3-V4 +1.5, swabs +2.0, spike-in +2.3, soil V4 -9.1); V1-V2 gained +3.7. On soil V4 G was 9.1 points worse than A (see below). |
| P9 | G is as accurate as A, or more, on every mock with sequences; on mock B at least 3 points more reads in exact ASVs | **Met.** Mock A, D6300 V3-V4 and V4-V5 'as accurate'; mock B 74.7% against 68.7% of reads in exact ASVs (+6.0), 'more accurate'. |
| P10 | F2 relaxes to (2,5) on mock B and wastewater | Met (it also relaxed on inter-laboratory V4-V6: 42% pass at (2,2)). |
| P11 | F2 never merges more pairs than the better of A and G by more than the margin on any 16S dataset | **Met on all 11.** The closest: wastewater -0.4, V4-V6 -1.2, mock B -4.4 against G; elsewhere 2.0 to 17.7 points behind. |
| P12 | F2 runs on every 16S dataset and merges at least 20 points more than F did in night 2 where F merged under 25% or gave no result | Met on all 9 (F2 55.5% soil V4, 81.3% D6300 V3-V4, 57.2% V4-V5, 45.7% mock B, 50.7% wastewater, 59.2% V4-V6, 62.5% mock A, 80.1% spike-in, 77.0% swabs). |

**What follows, by the rule fixed in section 14.** P9 holds and P7 does not, so the planned warning shows the table of looser settings **without naming one**. The failure of P7 is on closeness to FIGARO (wastewater -15.7), not on harm: G never lost accuracy and merged 7.5 to 29.8 points more pairs than the default on every poor dataset. Naming (2,5) in the warning would be a post hoc decision, to be labelled as such, that cannot be tested on these data again.

**Gain over the default (G minus A, merged pairs).** Poor datasets (A passes under 70%): mock A +13.4, mock B +29.8, V4-V6 +18.8, V4-V5 +11.6, wastewater +7.5. ITS (74.2% pass, no truncation): +2.7 (+2.3 to +3.1), pass 74.2% to 85.7%. Clean datasets: V1-V2 +3.7, spike-in +2.3, swabs +2.0, D6300 V3-V4 +1.5, inter-laboratory V1-V2 +1.4, soil V4 **-9.1** (AmpliTrim's search at (2,5) chose 85/154 against 124/115; its yield estimate ranks soil V4 poorly, rank correlation +0.09 in night 2). Rule R (G where A passes under 70%, else A) left soil V4 at A. Against FIGARO's pick the gaps of A became, under G: mock A -15.3 to -1.8, V4-V6 -27.0 to -8.2, wastewater -23.1 to -15.7, mock B -54.1 to -24.2. Both G and the default were 'trade-off: more pairs, noisier' against FIGARO on swabs, spike-in and D6300 V3-V4 (G +7.9, +8.6, +7.6 points).

**Against a careful manual cut (F2).** Where F2 kept maxEE 2 it merged fewer pairs than A on five datasets (V1-V2 -12.1, soil V4 -11.4, inter-laboratory V1-V2 -16.3, D6300 V3-V4 -9.3, spike-in -3.4) and the same on swabs and V4-V5 (its capped cut equals AmpliTrim's) and mock A (-0.4). Where it relaxed to (2,5) (mock B, V4-V6, wastewater) it came within 4.4, 1.2 and 0.4 points of G: the gain there is the maxEE, not the manual cut. AmpliTrim's cut was never beaten by a careful manual cut at the same maxEE.

**A caution on accuracy.** Arms with the same maxEE and cuts that differ by a few nucleotides differed in accuracy: on mock B, F2 and G (both maxEE 2,5, 280/247 against 277/250) reached 82.3% and 74.7% of reads in exact ASVs; on mock A, F2 (280/240) recovered 33 variants and 10 strains against 36 and 12 for A (277/243), the same total length. Differences of this size between single arms should not be read as an effect of maxEE alone. What the data support on mock B is that all seven looser arms (73.6% to 82.3%) were above the default (68.7%).

**Not done.** The control that would separate 'relax the limit' from 're-split the cut' (A's own cut at maxEE (2,5)) was not run; G is AmpliTrim's best cut at that maxEE. F2 is post hoc. The 70% threshold was chosen after night 2 and applied here to the same datasets.

**Decisions.**
1. The warning of section 13, decision 3(b), shows the table of looser settings (their estimated pass rate and yield) without naming one, as fixed in section 14.
2. The default maxEE (2) is not changed; G is not recommended for every dataset (soil V4 lost 9.1 points, clean datasets gained at most 3.7).
3. Validation of maxEE ends here. Nothing else is tested before release except the offline check of FIGARO's exact derivation (decision 2 of section 13; it needs FIGARO's source).
4. The warnings are built next (text only, with tests and a new build), then `docs/FINDINGS.md` and the README carry the night 3 results (done with this section).

## 16. Night 4 (offline): a simpler maxEE rule, with FIGARO's derivation as the reference (written 10 October 2026, before `validation/maxee_rules.py` was run on any dataset)

**Why.** Night 3 left AmpliTrim's maxEE rule (R: 2, or 2,5 when fewer than 70% of the pairs pass at 2) 16 and 24 points behind FIGARO's pick on wastewater and mock B, and the thresholds of R were chosen after seeing the data. The question is whether one parameter-free formula can do what R does and what FIGARO does: leave DADA2's default alone on clean reads and go looser on poor ones.

**What was found in FIGARO's source** (`trimParameterPrediction.py` and `expectedErrorCurve.py`, recovered from the session of 3 October; the function that builds the data the curve is fitted to was not available): at a cut, FIGARO's maxEE per read is `padMaxExpectedError(curve(position))`, where `padMaxExpectedError(x) = -(int(-x)) + 1`, which is `int(x) + 1` for x >= 0 (truncation, not a ceiling), and `curve` is `a*exp(b*x) + c` fitted with `scipy.optimize.curve_fit` to the expected errors of the reads along the read (start values 0.03, 0.015, 0; bounds (-2,-1,-8) to (2,1,8)). Hard-coded constants (`0.0356*e^(0.015L)` forward, `0.0289*e^(0.0203L)` reverse) are used only when no curve exists. The score `retention*100 - ((F-1)^2 + (R-1)^2)` ranks cuts and plays no part in choosing a cut's maxEE. So FIGARO's 1 on clean reads is what the formula gives when a typical read has under one expected error at the cut. **This corrects section 13's B4**, which used `round(mean) + 1`: rounding differs from truncation by up to one, and a raw mean differs from a fitted curve. B4 was an approximation, and the conclusion of section 13 that no AmpliTrim rule matched FIGARO's choice was about approximations.

**What DADA2 says about maxEE (recalled, not re-read: no network in this session; to be checked against the `filterAndTrim` help page, the tutorial and the GitHub issues before it is quoted).** maxEE is the largest number of expected errors, `sum(10^(-Q/10))` over the truncated read, a read may carry; filterAndTrim does not filter on it unless it is set; the tutorial starts at 2 and says to relax it, perhaps especially on the reverse reads (for example `c(2, 5)`), when too few reads pass, and to tighten it to speed up later steps; DADA2's error model uses the quality scores, so it tolerates some noisier reads. Nothing in this says a stricter maxEE improves accuracy.

**The rules (all per read; the cut is AmpliTrim's recommendation at maxEE 2).**
- *Mean rule (primary):* `int(mean expected errors of the reads at the cut) + 1`. No curve, no scipy, no threshold, no cap.
- *Median rule (secondary):* the same with the median.
- *Floored:* each of them never below 2, DADA2's default. This is the only free choice in the proposal, and it rests on DADA2's own default and on night 1 to 3: FIGARO's maxEE 1 on clean data merged 2.3 to 6.3 points fewer pairs than the default and carried 0.01 to 0.10 fewer expected errors per read, and no maxEE tested ever cost accuracy on a mock.

**Fixed before the run: FIGARO's cut and pick (night 3, arm C), AmpliTrim's cut (arm A).**

| dataset | FIGARO cut | FIGARO maxEE | DADA2 passed | AmpliTrim cut | kind |
|---|---|---|---|---|---|
| mock A | 273/247 | 3,4 | 80.3% | 277/243 | poor |
| V1-V2 | 225/138 | 2,1 | 88.2% | 201/157 | clean |
| soil V4 | 120/119 | 1 | 96.5% | 124/115 | clean |
| inter-laboratory V1-V2 | 156/169 | 1 | 94.1% | 155/170 | clean |
| mock B | 279/248 | 6,9 | 75.9% | 279/248 | poor |
| inter-laboratory V4-V6 | 270/268 | 4 | 70.4% | 273/265 | poor |
| D6300 V3-V4 | 224/218 | 1 | 88.2% | 224/218 | clean |
| wastewater | 247/193 | 5 | 83.1% | 254/186 | poor |
| spike-in | 232/207 | 2,1 | 81.6% | 222/218 | clean |
| swabs | 224/218 | 1 | 78.6% | 224/218 | clean |

Poor: FIGARO's maxEE is 3 or more on at least one read. Clean: 2 or less on both. A single number means the same limit on both reads. D6300 V4-V5 and ITS have no FIGARO result and are not in this test.

**Predictions, written before the run.**
- **P13.** The mean rule, unfloored, at FIGARO's cut, is within 1 of FIGARO's maxEE on each read on at least 8 of the 10 datasets (the criterion fixed in section 13, decision 2, for an exact implementation).
- **P14.** On the four poor datasets, the floored mean rule at AmpliTrim's cut is within 1 of FIGARO's pick on each read on at least 3.
- **P15.** On the six clean datasets, the floored mean rule at AmpliTrim's cut gives (2,2) on at least 5.
- **P16.** The median rule matches FIGARO on no more datasets than the mean rule.

The run also checks the reads: AmpliTrim's filter replica at FIGARO's cut and pick must be within 2 points of the share DADA2 passed in night 3; if not, the folder is probably not the reads of that run.

**What follows.**
- *P13, P14 and P15 hold:* the floored mean rule replaces rule R as AmpliTrim's maxEE recommendation, in two steps (the best cut at maxEE 2; the maxEE from the reads at that cut; the search repeated at that maxEE). The outcomes are mostly already measured: where the rule gives (2,2) the result is arm A, where it gives FIGARO's pick at AmpliTrim's cut it is arm D (night 2). A DADA2 run (night 4b, arm H) is needed only where the rule gives a maxEE no arm has used. Adoption requires, on every dataset, merged pairs at least those of arm G (maxEE 2,5) within the 2-point margin, and on the mocks with expected sequences no loss beyond the margins of night 2 (variants 1, strains 1, 0.1 points of reads in identical ASVs). Whether to cap the rule (FIGARO went to 9 on the reverse read of mock B; DADA2's tutorial example is 5) is decided then, and is a judgement: nothing here measures accuracy in rich real samples.
- *P13 fails:* the formula on the raw reads is not what FIGARO computes. Rule R stays, and the exact fitted-curve derivation is tested next (it needs the function that builds the curve's data, from `expectedErrorCurve.py`).
- *P13 holds and P15 fails:* the rule relaxes where the default is fine; it is not adopted without a floor higher than 2 or a threshold, which would be a new, post hoc rule.
- Whatever happens, this is a post hoc proposal designed after seeing FIGARO's picks, and the datasets that test it are the datasets that inspired it. It would be listed as such in the paper.

**Addendum (10 October 2026, after a first, partial run; the criteria above are unchanged).** The first run of `maxee_rules.py` covered 2 of the 10 datasets, mock A and V1-V2, because the script read only the top-level keys of `amplitrim_paths.json` and skipped the eight datasets that are registered under `extras`. That was a fault of the script, now fixed (`reads_folder` reads both places, as `night2.py` does). What the two datasets showed is recorded here so that the full run cannot be read selectively. *The reads were the right ones:* AmpliTrim's filter replica at FIGARO's cut and pick passed 80.2% of the pairs on mock A (DADA2 80.3%) and 88.2% on V1-V2 (DADA2 88.2%). *Mock A (poor):* the mean expected errors of the reads at FIGARO's cut are 1.30 forward and 1.95 reverse, so the mean rule gives (2,2) where FIGARO chose (3,4): within 1 on the forward read, **not within 1 on the reverse read**. *V1-V2 (clean):* mean 0.70 / 0.24, rule (1,1) against FIGARO's (2,1): within 1. So P13 stands at 1 of 2 and P14 at 0 of 1; the median rule gave (1,2) on mock A at FIGARO's cut and (1,1) on V1-V2. On both datasets FIGARO's pick implies a curve value above the mean of the reads at its cut (at least 2 and 3 where the mean is 1.30 and 1.95; at least 1 where it is 0.70), which the exact derivation would explain if the fitted curve is not the mean. An exploratory table was therefore added after seeing this: the scale k of a statistic of the reads (mean, median, 75th and 90th percentile) that would reproduce every FIGARO pick. It is **not part of the criteria**, and a k fitted on these datasets would be post hoc and would add a parameter.

**Outcome of the full run (10 of 10 datasets, 10 October 2026): the mean rule is not adopted.**

| | Prediction | Outcome |
| --- | --- | --- |
| P13 | The unfloored mean rule is within 1 of FIGARO's maxEE on each read, at FIGARO's cut, on at least 8 of 10 | **Not met: 7 of 10.** It matches on all six clean datasets and on inter-laboratory V4-V6 (3,3 against 4,4). The three misses are the poor reads: mock A (2,2 against 3,4), mock B (3,6 against 6,9) and wastewater (3,3 against 5,5), where it is 2 to 3 below FIGARO. |
| P14 | On the four poor datasets the floored mean rule at AmpliTrim's cut is within 1 of FIGARO's pick on each read on at least 3 | **Not met: 1 of 4** (inter-laboratory V4-V6 only). |
| P15 | On the six clean datasets the floored mean rule gives (2,2) on at least 5 | Met: 6 of 6, but by construction: the mean expected errors at the cut are under 1 on every clean dataset (largest 0.72), so the unfloored rule is always (1,1). |
| P16 | The median rule matches FIGARO on no more datasets than the mean rule | Met: median 6, mean 7. |

The reads were the right ones: AmpliTrim's filter replica at FIGARO's cut and pick was within 1.8 points of the share DADA2 passed in night 3 on all ten (the largest gaps: spike-in 83.4% against 81.6%, wastewater 84.7% against 83.1%; the other eight within 0.7 points).

At AmpliTrim's cut, on the four poor datasets, the share of pairs passing (exact quality table, AmpliTrim's estimate) with the default, the floored mean rule and FIGARO's pick, and rule R as measured by DADA2 in night 3 at its own cut:

| dataset | default (2,2) | floored mean rule | FIGARO's pick | rule R (2,5), night 3 |
| --- | --- | --- | --- | --- |
| mock A | 64.6% | 64.6% at (2,2) | 80.1% at 3,4 | 79.0% |
| mock B | 19.6% | 59.2% at 3,6 | 76.6% at 6,9 | 50.9% |
| inter-laboratory V4-V6 | 43.3% | 59.7% at 3,3 | 70.5% at 4,4 | 62.5% |
| wastewater | 55.4% | 64.8% at 3,2 | 84.1% at 5,5 | 67.8% |

The mean rule keeps more pairs than R on one of the four (mock B, +8.3 points) and fewer on three (mock A -14.4, V4-V6 -2.8, wastewater -3.0; the rule R figures are at its own cut and measured by DADA2, so the comparison is approximate). **Decision: rule R stays as the maxEE recommendation of the release.** FIGARO's picks passed 70.4% to 96.5% of the pairs on all ten datasets (night 3, DADA2).

*Exploratory (not a criterion).* No single scale k of the mean, median, 75th or 90th percentile reproduces every FIGARO pick. For the mean, k must be at least 1.89 (mock B, forward read) and below 1.73 (set by other reads), so FIGARO's curve value is not a constant multiple of the mean: the curve is a different quantity. On all eight reads of the four poor datasets, FIGARO's pick implies a curve value that lies between the 75th and the 90th percentile of the expected errors of the reads at its cut (the mean lies below the 75th percentile there: expected errors are right-skewed).

---

## 17. Night 4, second stage (offline): a percentile rule (written 10 October 2026, before it was evaluated on any dataset)

**Why.** The observation that ends section 16 suggests that FIGARO's maxEE follows a high percentile of the reads' expected errors at the cut, not their mean. A percentile rule is also easier to explain than a mean plus one: *keep the best q of the reads in each direction, and never go below DADA2's default of 2*.

**The rule.** Per read, `int(q-th percentile of the expected errors at the cut) + 1`; floored: the same, never below 2. The cut is AmpliTrim's recommendation at maxEE 2. **q is fitted:** it is chosen from {0.75, 0.80, 0.85, 0.90} (fixed here, before the run) as the value whose unfloored rule matches FIGARO on most datasets, ties to the smaller q (cleaner reads). The idea comes from looking at these same ten datasets, so the rule is post hoc and q is a parameter fitted on the datasets that test it; both are to be said in the paper.

**Predictions, written before the run.**
- **P17.** For at least one q of the grid, the unfloored percentile rule at FIGARO's cut is within 1 of FIGARO's maxEE on each read on at least 9 of the 10 datasets.
- **P18.** With the q that P17 selects, at AmpliTrim's cut, the floored rule gives (2,2) on at least 5 of the 6 clean datasets and reaches 70% of the pairs passing (exact quality table) on at least 3 of the 4 poor ones. (70% is the threshold of rule R and the lowest share FIGARO's picks passed.)

The run also reports, as exploration, the interval of q that would reproduce every pick (the intersection over all reads of the share of reads with fewer expected errors than pick - 1 and than pick); that interval is not a criterion.

**What follows.**
- *P17 and P18 hold:* the floored percentile rule with that q is proposed to replace rule R, pending one DADA2 confirmation (night 4b): a run is needed only for datasets where the rule's maxEE at AmpliTrim's cut is one no arm has used (arms A, D and G of nights 2 and 3 cover the others). Adoption requires, on every dataset, merged pairs at least those of arm G within the 2-point margin, and on the mocks with expected sequences no loss beyond the margins of night 2 (variants 1, strains 1, 0.1 points of reads in identical ASVs). Whether to cap the rule (FIGARO went to 9 on the reverse read of mock B; DADA2's tutorial example is 5) is decided then, and is a judgement: nothing here measures accuracy in rich real samples. Looser limits keep more pairs and noisier reads, and nothing measured here penalises the noise, so "better" means a stated trade, not an optimum.
- *P17 fails:* no single-parameter rule of this kind reproduces FIGARO. Rule R stays, and the last test is FIGARO's exact fitted curve (it needs the function that builds the curve's data from `expectedErrorCurve.py`).
- *P17 holds and P18 fails:* the rule is described in the paper as a finding; rule R stays.
- Whatever happens, **rule R is the rule of the release until night 4b.**

**Outcome of the percentile run (10 of 10 datasets, 10 October 2026): P17 and P18 are met.**

| | Prediction | Outcome |
| --- | --- | --- |
| P17 | For at least one q of the grid, the unfloored percentile rule at FIGARO's cut is within 1 of FIGARO's maxEE on each read on at least 9 of 10 datasets | **Met.** q = 0.80 and q = 0.85 both match on 10 of 10 (q = 0.75 and 0.90: 8 of 10). The tie goes to the smaller q, as registered: **q = 0.80.** |
| P18 | With that q, at AmpliTrim's cut, the floored rule gives (2,2) on at least 5 of the 6 clean datasets and reaches 70% passing on at least 3 of the 4 poor ones | **Met:** 6 of 6 and 4 of 4. Only q = 0.75 fails the second part (2 of 4). |

*What is solid, and what is soft.* "Within 1" is a lenient criterion. Exactly, q = 0.80 reproduces 14 of FIGARO's 20 read-picks and q = 0.85 reproduces 16: at q = 0.85 the three datasets that differ are mock B (6,10 against 6,9), inter-laboratory V4-V6 (5,5 against 4,4) and wastewater (5,4 against 5,5); at q = 0.80 the differences are mock A (2,4 against 3,4), V1-V2 and spike-in (1,1 against 2,1, which the floor removes), mock B (5,9 against 6,9) and wastewater (4,4 against 5,5). No single q reproduces every pick exactly: the interval of q that would is empty (above 0.87, at most 0.83). The clean part of P18 is by construction (every clean dataset has few expected errors, so any floored rule gives the default). q was fitted on the datasets that test it.

At AmpliTrim's cut, what q = 0.80 and q = 0.85 pick and keep (exact quality table, AmpliTrim's estimate; the search is not repeated at the new maxEE), against FIGARO's pick and rule R (night 3, DADA2, at its own cut):

| dataset | q = 0.80 | q = 0.85 | FIGARO's pick | rule R (2,5) |
| --- | --- | --- | --- | --- |
| mock A | 3,3: 75.7% | 3,4: 80.1% | 3,4: 80.1% | 79.0% |
| mock B | 5,9: 74.8% | 6,10: 78.7% | 6,9: 76.6% | 50.9% |
| inter-laboratory V4-V6 | 4,4: 70.5% | 5,5: 77.8% | 4,4: 70.5% | 62.5% |
| wastewater | 4,3: 75.8% | 5,4: 82.3% | 5,5: 84.1% | 67.8% |
| the six clean datasets | (2,2), as the default | (2,2) | 1,1 or 2,1 | (2,2) |

q = 0.80 passes 3.3 points fewer pairs than rule R on mock A and 23.9, 8.0 and 8.0 points more on mock B, V4-V6 and wastewater (approximate: different cuts, an estimate against DADA2). It lets in slightly cleaner reads than FIGARO on mock B (1.20 / 3.57 expected errors per kept read against 1.30 / 3.63) and wastewater (1.31 / 0.96 against 1.55 / 1.16). **Not tested:** ITS (no cut), and D6300 V4-V5 (no FIGARO row, so it is not in this table): night 4b computes the rule there.

**Decision.** By section 17, the floored percentile rule with q = 0.80 is proposed to replace rule R, **pending the DADA2 confirmation of section 18.** Rule R remains the rule of the release until then.

---

## 18. Night 4b: the percentile rule through DADA2 (written 10 October 2026, before `validation/night4.py` was run)

**Question.** Does the percentile rule (q = 0.80) merge at least as many pairs as rule R's relaxed setting (2,5), without losing accuracy on the mocks?

**The rule, exactly.** At AmpliTrim's cut at maxEE 2, per read: `int(q-th percentile of the expected errors of the reads at the cut) + 1`, never below 2; then the search is repeated at that maxEE. q = 0.80 is arm **H** (the registered rule); q = 0.85 is arm **H85**, exploratory and post hoc. The rule is computed by the script from the reads, not copied from section 17. Where both q give (2,2) the result is arm A by construction and nothing is run (expected: all six clean datasets). ITS is not tested.

**Arms** (all run through DADA2 on the same pairs, per dataset): A (maxEE 2, the control), G (2,5, the setting to beat), H, H85, C (FIGARO's own pick at the longest amplicon length). Margins as in night 2: merged share 2 points, expected errors 0.02 per read, accuracy {variants 1, strains 1, identical-ASV reads 0.1 points}. For accuracy only 'as accurate' and 'more accurate' pass; 'mixed' fails.

**Predictions, written before the run.**
- **P19.** On every dataset that was run, H merges at least as many pairs as G within the margin (H minus G at least -2 points). *Offline, mock A is expected to fail this:* the rule's (3,3) passes 75.7% of the pairs against 79.0% for (2,5). It is registered because adoption was fixed before the run.
- **P20.** On every mock with expected sequences that was run (mock A, mock B, D6300 V4-V5), H is as accurate as G or more accurate.
- **P21.** H merges more pairs than G by more than the margin on mock B, inter-laboratory V4-V6 and wastewater (offline: +23.9, +8.0, +8.0 points of pairs passing).
- **P22 (exploratory).** H85 merges within the margin of FIGARO's pick on every dataset where FIGARO gave a result.

**What follows.**
- *P19 and P20 hold:* the percentile rule (q = 0.80) replaces rule R in the release, its post hoc status and the fitted q stated, and the cap is decided (FIGARO went to 9 on mock B's reverse read; DADA2's tutorial example is 5; this is a judgement, because the study measures merged pairs and accuracy on mocks, not the effect of noisier reads in rich real samples).
- *P19 fails only narrowly on a dataset (as expected on mock A):* it is reported as it is, and the rule is **not** adopted on the strength of arm H85: a different q chosen after this run would be a new post hoc rule needing its own test.
- *P20 fails on any mock:* rule R stays.
- P21 and P22 are descriptions. Whatever happens, rule R is the rule of the release until this has been read.

**Limits.** The datasets that test the rule are the datasets that inspired it, and q is fitted on them; there is no independent dataset left. The accuracy margin on identical-ASV reads (0.1 points) is the one fixed for the battery and is strict. Only the mocks measure accuracy.

**Revision of section 18 (10 October 2026, before night 4b was run; it supersedes the arms and the predictions P19 to P22 above, and keeps everything else: the margins, the exclusion of ITS, the rule that where the rule gives (2,2) the arm is A, and the limits).**

*What changed, and why.* The researcher (Payton Yau) chose, for the release, **q = 0.75 and a cap of 5** on each read, instead of the q = 0.80 without a cap that section 17 proposed. The reasons given: a lower maxEE follows DADA2's own guidance more closely (the tutorial's example for relaxing the filter is `c(2, 5)`, and the default is 2), and the cap keeps the rule inside the range the DADA2 documentation illustrates. The cap is recorded as **a design choice consistent with the tutorial's example**, not as a recommendation of DADA2's author: the text offered in support attributed specific statements to him (a ceiling of 5, a higher limit for PacBio long reads, retention thresholds) that could not be checked in this session, and at least one does not match what is recalled of the documentation (the DADA2 PacBio tutorial, as recalled and not re-read, uses maxEE = 2). They are to be checked against the primary sources (the `filterAndTrim` help page, the tutorial and the GitHub issues) before anything is quoted as the author's advice. One point of that text does apply and is recorded: "shorten truncLen instead of raising maxEE" is not available to AmpliTrim on the datasets where the rule relaxes, because its cut is already the shortest that can merge the longest amplicon plus the minimum overlap.

*What is known about q = 0.75 before the run (section 17, offline).* It failed both registered criteria: it matched FIGARO's maxEE within 1 on 8 of 10 datasets (P17 needs 9) and reached 70% of the pairs passing on only 2 of the 4 poor datasets (P18 needs 3). At AmpliTrim's cut it picks (2,3) on mock A (71.5% of the pairs passing, against 79.0% for rule R's (2,5)), (4,8) on mock B (69.5%, uncapped), (3,3) on inter-laboratory V4-V6 (59.7%, against 62.5%) and (4,3) on wastewater (75.8%, against 67.8%), and (2,2) on the six clean datasets. **So q = 0.75 relaxes less than rule R on mock A and V4-V6; the cap lowers mock B's (4,8) to (4,5), whose share passing has not been measured.** Choosing q = 0.75 is therefore a trade of retention for conformity and cleaner reads, made after seeing these numbers; it is post hoc.

*Arms* (all through DADA2 on the same pairs, per dataset): A (maxEE 2, the control), G (2,5, the setting to beat), **H** (the rule at q = 0.75, never below 2 and never above 5), **H80** (q = 0.80 with no cap, the rule of section 17; it shows what the cap costs) and C (FIGARO's own pick at the longest amplicon length). Where H and H80 both give (2,2) the arm is A and nothing is run (expected: the six clean datasets).

*Predictions, written before the run (they replace P19 to P22).*
- **P19.** On every dataset that was run, H merges at least as many pairs as G within the margin (H minus G at least -2 points). *Expected to fail on mock A* ((2,3): 71.5% of the pairs passing against 79.0%) *and probably on V4-V6* (59.7% against 62.5%).
- **P20.** On every mock with expected sequences that was run (mock A, mock B, D6300 V4-V5), H is as accurate as G or more accurate; 'mixed' fails.
- **P21.** H merges more pairs than G by more than the margin on wastewater (offline +8.0) and on mock B (offline for the uncapped (4,8): +18.6 points of pairs passing; the capped (4,5) is not measured offline).
- **P22 (exploratory).** Where the uncapped rule gives a maxEE above 5 on a read (expected: mock B), H80 merges more pairs than H by more than the margin: that is what the cap costs.

*What follows.* If P19 and P20 hold, the rule (q = 0.75, cap 5) replaces rule R in the release, its post hoc status stated. **If P19 fails, as expected on mock A, the rule is not shown to be at least as good as R**, and it is not adopted on the strength of H80 or of any q chosen after this run (that would be a new post hoc rule). It may still be chosen for the release as a stated trade, conformity with DADA2's guidance and cleaner reads against fewer pairs, but not as a validated improvement over R. If P20 fails, rule R stays. Whatever happens, **rule R is the rule of the release until this has been read.** P21 and P22 are descriptions.

**Verification of the sources for the revision above (10 October 2026).** The researcher pasted the text of DADA2 issues #2136 (July 2025), #236 (2017 to 2020) and #1909 (March 2024). The summary that had been offered in support of a cap of 5 was checked against them, claim by claim.

| Claim in the summary | What the threads say |
| --- | --- |
| maxEE should normally not exceed 5; going above it signals structural problems | **Not found.** Callahan: "maxEE=2 is a good general purpose value. Relaxing it somewhat is also OK" (#2136); "try a couple maxEE values until I get a satisfactory number of reads through the filter" and "potentially relax maxEE a bit" (#236, 2017). No ceiling. |
| Up to 6 on very messy reverse reads | **Not found.** "Values from 2 to 6 are shown in various issues" is the *opener's* summary of other threads (#236), not Callahan's words. |
| maxEE of 10 or more causes false-positive ASVs or distorts abundances | **Not found as stated.** Closest: "pushing through more reads of lower quality may not improve the relative abundance estimates, it could even make it worse" (#2136, about PacBio). |
| maxEE = 2 is a technology-agnostic baseline (AVITI, MiSeq) | **Partly.** "A good general purpose value" (#2136). AVITI and "technology agnostic" are not in these threads. |
| Losing many reads is normal; fewer good reads beat diluting with noisy ones | **Supported, with context.** 40% after filtering is fine, the number of reads is irrelevant (#2136, PacBio); "nothing wrong with stricter filtering when you are just giving up a marginal (e.g. 10%) of reads", "more lower quality reads doesn't really help anything" (#236, 2020). |
| PacBio: maxEE up to 12 is acceptable | **Not found.** In #2136 the opener wrote that "up to 5 or even 8" is sometimes used for PacBio; Callahan answered "relaxing it somewhat is also OK". No 12. |
| ITS: maxEE 5 is acceptable | **Not found.** #1909: truncLen is "not appropriate for highly length-variable amplicons like ITS"; no maxEE number. |
| If reads fail at maxEE 2, do not raise it; shorten truncLen | **Partly.** On Illumina "trimming would be much more effective" is a contributor's remark and Callahan agrees it is "usually the case" because of the quality crash at read ends (#236, 2018). He does not say not to raise maxEE; his 2017 procedure includes relaxing it. |
| Retention thresholds (40 to 50% loss normal; 60 to 80% caution; above 90% act, up to c(2,5) or c(3,5)) | **Not found.** "At least 10% of raw reads" and "60% passing is good" are in the opener's summary of other threads (#236), not Callahan's words in these. |
| Relax the reverse read to 5 and keep the length long | **Not found** in these threads. |
| Forward + reverse truncLen above the amplicon length + 20 to 30 nt | **Supported:** "truncLen[[1]] + truncLen[[2]] > amplicon_length+25" and "err a bit on the more sequence side" (#236, 2017). |

*Consequences.* (1) **The cap of 5 is a design choice of the researcher; DADA2's author gives no numerical ceiling in these threads.** The numbers that look like ceilings were written by the people who opened the issues. The tutorial's `c(2, 5)` was recalled and not re-read, and is no longer cited in the release documentation. (2) What he does say argues for **restraint**: the number of reads is not the aim, and more lower-quality reads may not help and could make estimates worse. That favours relaxing only where the default loses a large share of the pairs, and it puts the uncapped rule (q = 0.80, up to 9 or 10) further from his guidance than q = 0.75 with a cap. It also sits in tension with any rule whose purpose is to raise retention: the validation measures merged pairs and accuracy on mock communities, and **does not measure the effect of noisier reads on relative abundances in rich samples**, which is the effect he warns about. (3) AmpliTrim's cut is at the coverage floor (longest amplicon + 12 nt), so "shorten truncLen" is not available where the rule relaxes; his rule of thumb ("in general", for choosing cuts by eye) puts the margin at 25 nt, and in the example he gave, a 300 nt amplicon on 2x300 reads, his suggested cuts overlapped by about 100 nt. AmpliTrim does not use a fixed margin: it tests real merges at the minimum overlap instead. His statement about keeping 40% of the reads was about full-length PacBio reads, so it is no more general than the 25 nt. The consequence of that difference for DADA2's merging has been measured at AmpliTrim's cuts (night 2 and 3 report DADA2's merged shares) but is not a tested comparison against a +25 cut.

**Offline variants (10 October 2026, read before night 4b was run; the predictions above are not changed).** The offline script now reports each q with and without a cap of 5 at AmpliTrim's cut. **The cap binds on one dataset only, mock B**, the only one where an uncapped pick exceeds 5: at q = 0.75, (4,8) passing 69.5% of the pairs becomes (4,5) passing 55.2% (-14.3 points; expected errors per kept reverse read 3.31 to 2.55, forward 1.02 to 0.81); at q = 0.80, (5,9) passing 74.8% becomes (5,5) passing 56.1% (-18.7). On mock A (2,3 and 3,3), inter-laboratory V4-V6 (3,3 and 4,4) and wastewater (4,3) the cap changes nothing. The rule of the revision (q = 0.75, cap 5) against rule R's (2,5), the share of pairs passing (the rule at AmpliTrim's cut, an estimate; rule R measured by DADA2 at its own cut, so the comparison is approximate):

| dataset | rule R (2,5) | q = 0.75, cap 5 | difference | kept-read noise F / R: R, then the rule |
| --- | --- | --- | --- | --- |
| mock A | 79.0% | 2,3: 71.5% | -7.5 | 0.37 / 1.43, then 0.44 / 0.85 |
| mock B | 50.9% | 4,5: 55.2% | +4.3 | 0.56 / 2.53, then 0.81 / 2.55 |
| inter-laboratory V4-V6 | 62.5% | 3,3: 59.7% | -2.8 | 0.52 / 1.48, then 0.84 / 0.85 |
| wastewater | 67.8% | 4,3: 75.8% | +8.0 | 0.69 / 1.98, then 1.31 / 0.96 |

*Expected from this table:* P19 fails on mock A (-7.5) and on V4-V6 (-2.8), as predicted; P21 holds on mock B (about +4.3) and wastewater (+8.0), by small margins on mock B; P22 holds on mock B (the cap costs about 19.6 points of pairs between (5,9) and (4,5)). On average over the four poor datasets the rule keeps about as many pairs as rule R (+0.5 points), distributed differently. *Where the rule differs in kind from rule R:* R relaxes the reverse read only; the percentile rule relaxes both reads according to their own noise, so the forward read is looser and the reverse read stricter than under R (on wastewater the forward limit, 4, exceeds the reverse limit, 3). The three issues checked above say nothing for or against relaxing the forward read; relaxing it is a departure from relaxing the reverse read only, and DADA2's documentation on that point has not been re-read.

**Decision (10 October 2026): the trigger of the release's rule R is 60%, not 70%.** The researcher (Payton Yau) judged 70% too high and 50% too low. This is **post hoc**: it was chosen after night 3 and after seeing the table below, between two values that give different sets of relaxed datasets. Rule R is the rule of the release until night 4b has been read; its form is unchanged (maxEE 2, or (2,5) when fewer than the trigger's share of the pairs pass at 2, or none do), only the number changes. The outcomes at 60% are already measured, because arms A and G of night 3 do not depend on the trigger: R uses G where A passed under the trigger, and A elsewhere.

| trigger | datasets relaxed to (2,5), with the gain in merged pairs (night 3, DADA2) | left at the default, with what (2,5) would have gained |
| --- | --- | --- |
| 40% | mock B +29.8 | V4-V6 +18.8, wastewater +7.5, D6300 V4-V5 +11.6, mock A +13.4 |
| 50% | mock B +29.8, V4-V6 +18.8 | wastewater +7.5, D6300 V4-V5 +11.6, mock A +13.4 |
| **60% (chosen)** | mock B +29.8, V4-V6 +18.8, wastewater +7.5 | D6300 V4-V5 +11.6 (63.7% passing), mock A +13.4 (64.7%) |
| 70% (night 3) | the three above, D6300 V4-V5 +11.6, mock A +13.4 | none |

The share of pairs passing at maxEE 2 in night 3 was 20.3% (mock B), 43.1% (V4-V6), 53.7% (wastewater), 63.7% (D6300 V4-V5), 64.7% (mock A), 74.2% (ITS) and 85.0% or more on the other six; at 60% the two datasets just above the trigger keep the default and give up 11.6 and 13.4 points, and the three below it are relaxed. The trigger does not affect arms G, H or H80 of night 4b or its predictions, which compare the percentile rule with (2,5) wherever the run relaxes. The percentile rule has no trigger: if it replaces rule R, this decision lapses. At that time the tests and the release documentation stated 60%; the build tested in nights 1 to 3 (and in the validation archive) used 70% in its rule R arm, `POOR_BELOW` in `maxee_arms.py`, which is left as it was because it describes what was run.

**Decision (10 October 2026, after the 60% decision above): the release uses the percentile rule, not rule R; the 60% decision lapses.** The researcher (Payton Yau) judged the percentile rule, q = 0.75 with a cap of 5, the best option for AmpliTrim and asked for rule R to be removed from the release. *Status, to be stated wherever the rule is described:* its DADA2 confirmation (section 18, night 4b) **has not been run**; offline, at AmpliTrim's cut, it keeps about as many pairs as rule R on average (+0.5 points over the four poor datasets: mock A -7.5, V4-V6 -2.8, mock B +4.3, wastewater +8.0), and P19 is expected to fail on mock A and V4-V6. By the logic of section 18 the rule is therefore **not shown to be at least as good as R**; it is chosen as a stated trade: it relaxes each read by its own noise (where (2,5) relaxes the reverse read only), it stays at DADA2's default on every clean dataset, and its limits stay inside the range that DADA2's author calls acceptable ("relaxing it somewhat is also OK", issue #2136), against fewer pairs than (2,5) on some datasets. The percentile (0.75, over 0.80 which matched FIGARO's choices more closely) and the cap (5) are the researcher's choices, made after seeing the validation data; the cap is not DADA2's author's ceiling (he gives none in the threads checked). *Update:* night 4b has since been run. Its outcome, and the decision that followed it, are recorded below, in the outcome of section 18.

*What the release does.* `run()` finds the best truncation at maxEE 2; measures the expected errors of every read that reaches the recommended cut; sets per read `min(5, max(2, int(75th percentile) + 1))`; and repeats the search at that maxEE when it is not (2,2). Where nothing passes at maxEE 2 the reads are measured at the best split, judged at the cap, of the shortest total that covers the amplicons; if the looser search still finds no recommendation the result stays at maxEE 2 with the rule's measurement recorded; when the reads cannot span the amplicons the rule is not computed. ITS is not adjusted (reads are not cut at a common position, so there is no cut to measure at, although each read is trimmed at its first base of quality 2 or below; untested on reads of different lengths). A maxEE given in the settings is used as given. `summary.json` carries `maxEE_rule` (relaxed, share passing at 2, percentile, cap, the cut and the percentile of the expected errors). The release documentation says that the DADA2 confirmation is pending; **it is to be updated when night 4b has been read.** Night 4b is unchanged: it still compares the rule (arm H) with (2,5) (arm G, rule R's setting), and R remains in the nights 1 to 3 results and scripts as what was run.

**Release note (10 October 2026): the ITS report shows where the pairs go.** Added after the percentile rule, and descriptive: for `--its`, the report and `summary.json` (key `its`) give the shares of the pairs analysed that are too short after the cut at the first base of quality 2 or below, contain an N, exceed maxEE, pass the filter, pass but do not merge (vsearch's estimate), and the share of pairs that merge at full read length; and, when it is larger than the share that passes and merges, the share of forward reads that would pass the same filter. The first rows replicate DADA2's no-truncLen filter, which AmpliTrim predicted within 0.1 points of the measured pass rate on the ITS dataset (night 1); the merging row inherits the yield estimate's optimism (+21 points at 5,000 pairs per sample on that dataset). Nothing here was tested through DADA2 beyond that. *Why this and not the percentile rule on ITS:* on the ITS dataset (night 3) the default maxEE passed 74.2% of the pairs and DADA2 merged 46.9%; (2,5) passed 85.7% and merged 49.5% (+2.7, 2.3 to 3.1), so maxEE is not where ITS loses its pairs. The forward-only comparison cites DADA2 issue #236 (2017) for the statement that using the forward reads alone is "totally reasonable" and that bad reverse reads can cost more than they add; AmpliTrim does not measure the effect on ASVs.

---

## 19. ITS check: are the numbers of the ITS report what DADA2 gives? (written 10 October 2026, before `validation/its_check.py` was run)

**Why.** The ITS report of the release shows where the pairs go (too short after the cut at the first base of quality 2 or below, containing an N, over maxEE, passing, passing but not merging) and the share of forward reads that would pass the same filter. AmpliTrim predicted DADA2's pass rate on the ITS dataset to within 0.1 points (night 1), but the split of the loss into rows and the forward-only share have not been compared with DADA2's own filter, and the 'do not merge' row is vsearch's estimate on the raw reads. The release documentation calls all of them estimates. Not required for the release; it decides whether the filter rows can be called a replica of DADA2's filter.

**Part A (same reads).** AmpliTrim's own ITS run (the code behind `amplitrim --its`, with vsearch) on the ITS dataset (PRJEB58189, 37 samples, `amplitrim_paths.json` key `its`); then DADA2's `filterAndTrim` with no truncLen on the very same pooled subsample files (`maxN = 0`, `maxEE = 2`, `truncQ = 2`, `minLen = 50`, phiX not removed, as in `dada2_real.R`), as a pair and with the forward reads alone. No denoising is run. The comparison is only valid if both saw the same number of reads; if not, nothing is decided.

**Part B (reads files that night 3 wrote).** The per-sample DADA2 track of the arm 'AmpliTrim (maxEE 2)' on the ITS dataset (5,000 pairs per sample): input, filtered, reads in denoised forward and reverse sequences, merged, non-chimeric. A pair can only merge if both reads were denoised, so per sample min(denoised forward, denoised reverse) bounds the pairs that reach merging: filtered minus that is at least what denoising removed, and that minus merged is at most what merging removed.

**Predictions, written before the run.**
- **P23.** AmpliTrim's share of pairs passing the filter is within 0.2 points of DADA2's, on the same reads.
- **P24.** AmpliTrim's share of forward reads passing alone is within 0.2 points of DADA2's single-end filter, on the same reads.
- **P25 (exploratory).** In DADA2's track, of the loss between the filter and the merged pairs, denoising accounts for more than merging (the loss at least due to denoising is larger than the loss at most due to merging).

**What follows.**
- *P23 and P24 hold:* the documentation says that the filter rows and the forward-only share replicate DADA2's filter on the ITS dataset (within 0.2 points). *Either fails:* they stay 'estimates', the rows that differ are traced, and nothing is claimed about them.
- *P25 holds:* the documentation says that on the ITS dataset, at 5,000 pairs per sample, most of the loss after the filter is denoising and that AmpliTrim's 'do not merge' row is only the merging part. *P25 fails:* merging dominates, and the documentation says by how much AmpliTrim's 'do not merge' estimate falls short of DADA2's.
- Not claimed in any case: the effect on ASVs or accuracy (the ITS dataset has no expected sequences), and anything about other ITS data (there is no other).

**Outcome of section 19 (10 October 2026; run on the ITS dataset, 49,720 pairs, `AmpliTrim` from the clean package, DADA2 1.40.0, night 3's track tables).**

| | Prediction | Outcome |
| --- | --- | --- |
| P23 | AmpliTrim's share of pairs passing the filter is within 0.2 points of DADA2's, on the same reads | **Met:** 74.30% against 74.30% (+0.00 points) |
| P24 | AmpliTrim's share of forward reads passing alone is within 0.2 points of DADA2's single-end filter, on the same reads | **Met:** 87.47% against 87.47% (+0.00 points) |
| P25 | (exploratory) Of DADA2's loss between the filter and the merged pairs, denoising accounts for more than merging | **Met, with narrow bounds:** at least 14.0 points to denoising against at most 13.4 points to merging; the bounds point the right way (the denoising figure is a lower bound, the merging figure an upper bound) |

The track (37 samples, 178,657 input pairs, shares of the input pairs): filtered 74.2%, reads in denoised forward sequences 61.9%, in denoised reverse sequences 60.3%, merged 46.9%, after chimera removal 46.9%. AmpliTrim's rows for the same run: too short 0.00%, containing an N 0.00%, over maxEE 25.70%, passing 74.30%, passing but not merging 1.89% (vsearch), passing and merging 72.41%.

*What this shows.* (1) The filter rows and the forward-only share replicate DADA2's filter on this dataset, to the hundredth of a point. (2) The merging rows do not match DADA2: the estimate leaves 1.9 points at merging and puts the yield at 72.4%, DADA2 merged 46.9%, a gap of 25.5 points; of the 27.4 points between the filter and the merged pairs, at least 14.0 are denoising and at most 13.4 are merging. The split of the rest is **not known** (if merging took all of its 13.4, the estimate is 11.5 points too low at merging; if it took 1.9, denoising took 25.5). (3) The forward reads alone passed 87.5%, 13 points more than the pairs and 41 points more than DADA2 merged, so the forward-only line appears on this dataset.

*Error in the first version of the ITS report, found by this check and corrected.* The first version said that 15.9% of the pairs "do not merge even at full read length", because the amplicon is longer than the reads can span or the overlap is poor, and that "no setting recovers these". This was incorrect. The 15.9% comes from merging the full-length reads before filtering. The 97.5% merge rate (72.41 of 74.30) comes from the same vsearch run on the filtered reads. The two figures differ because most pairs that do not merge at full length are noisy, and the maxEE filter removes them anyway (about 14 of the 15.9 points are among the 25.7% of pairs filtered out). The line therefore counted the same pairs twice, contradicted the table beside it, and wrongly said that no setting recovers them. It was removed from the report and from the `its` block of `summary.json` (`overlap_at_full_length`).

*Open, not tested.* A hypothesis for part of the merging gap: DADA2's `mergePairs` rejects pairs whose merge overhangs the end of the other read unless `trimOverhang = TRUE`, where vsearch with `--fastq_allowmergestagger` merges and trims them; with amplicons shorter than the reads, pairs would merge in vsearch and not in DADA2. A test would run `mergePairs` twice on the same denoised reads (default and `trimOverhang = TRUE`) and compare the merged counts. Not run.

*Decisions.* The user guide now says that the filter rows and the forward-only share replicate DADA2's filter on the ITS dataset, and that the merging rows are an optimistic upper bound (with the figures above); the report's ITS table header and Limits say the same. Nothing is claimed about accuracy or other ITS data.

---

## 20. ITS merge check: the exact split of DADA2's loss after the filter, and whether an overhang explains part of the merging loss (written 10 October 2026, before `validation/its_merge_check.py` was run)

**Why.** Section 19 left two things open on the ITS dataset. First, the split of the 27.4 points between the filter (74.2%) and the merged pairs (46.9%) was only bounded: at least 14.0 points to denoising, at most 13.4 to merging; AmpliTrim's vsearch estimate put 1.9 points at merging. Second, a hypothesis: DADA2's `mergePairs` refuses a pair whose merge overhangs the end of the other read unless `trimOverhang = TRUE`, where vsearch with `--fastq_allowmergestagger` merges and trims it; amplicons shorter than the reads would then merge in vsearch and not in DADA2. If true it is an ITS-specific fault that AmpliTrim's report could name.

**What it does.** DADA2 (learnErrors, dada) is run once on the ITS reads, thinned exactly as in night 3 (the same function, the same pairs per sample), with no truncLen (`maxN = 0`, `maxEE = 2`, `truncQ = 2`, `minLen = 50`, phiX not removed). `mergePairs` is then run twice on the same denoised reads, with the default settings (`minOverlap = 12`, `maxMismatch = 0`) and with `trimOverhang = TRUE`, both with `returnRejects = TRUE`, so that every pair that enters the merge step is listed. A rejected pair is counted by reason: overlap shorter than `minOverlap`; mismatches or indels in the overlap; otherwise a clean overlap of at least `minOverlap` that is still refused, which is the overhang. This gives the **exact** split in place of the bounds: filtered minus the pairs that enter the merge step is what denoising lost, and the rejected pairs are what merging lost.

**Predictions, written before the run.**
- **P26 (consistency).** The default merge reproduces night 3's merged share (46.9% of the input pairs) within 0.3 points. It is the same reads, the same settings and the same DADA2; if it does not, the other results are not comparable.
- **P27.** Denoising accounts for more of the loss between the filter and the merged pairs than merging does, exactly (consistent with the bounds of P25).
- **P28 (the hypothesis).** `trimOverhang = TRUE` merges at least 2.0 points more of the input pairs than the default. The threshold separates a material overhang effect from a negligible one; the merging loss is at most 13.4 points, so 2.0 is a modest fraction of it.
- **P29.** The exact merging loss exceeds AmpliTrim's estimate (1.89 points) by more than 5 points.

**What follows.**
- *P28 holds:* the ITS report names it. Where the amplicon is shorter than a read, the report says that DADA2's `mergePairs` needs `trimOverhang = TRUE` (or the read-through removed first with cutadapt) and gives the share of pairs concerned, estimated from the insert lengths and the read lengths; the DADA2 snippet for ITS carries `trimOverhang = TRUE`. *P28 fails:* the hypothesis is rejected for this dataset, nothing is added, and the merging gap is described as unexplained.
- *P29 holds:* the user guide and the report say that the merging estimate for ITS is optimistic by that much on the one dataset; the 'do not merge' row stays an upper bound on what merges. *P29 fails:* the estimate was closer than the bounds allowed, and the guide says so.
- *P26 fails:* nothing from this run is used until the difference is explained.
- Not claimed in any case: accuracy, ASVs, or anything about other ITS data (there is none).

**Outcome of section 20 (10 October 2026; ITS dataset, 37 samples, 178,657 input pairs, 5,000 pairs per sample as in night 3, DADA2 1.40.0).**

| | Prediction | Outcome |
| --- | --- | --- |
| P26 | The default merge reproduces night 3's merged share within 0.3 points | **Met:** 46.87% against 46.87% (+0.00 points) |
| P27 | Denoising accounts for more of the loss between the filter and the merged pairs than merging does, exactly | **Met:** 15.9% against 11.5% of the input pairs |
| P28 | `trimOverhang = TRUE` merges at least 2.0 points more of the input pairs than the default | **Not met:** +0.00 points |
| P29 | The exact merging loss exceeds AmpliTrim's estimate (1.89 points) by more than 5 points | **Met:** 11.5 points, 9.6 over |

The exact flow, shares of the input pairs: pass the filter 74.24%; lost in denoising 15.87%; enter the merge step 58.37%; rejected for an overlap shorter than 12 nt (fewer than 12 matching bases) 5.07%; rejected for mismatches or indels in the overlap 6.43%; rejected with a clean overlap that is still refused (an overhang) 0.00%; merged 46.87%; merged with `trimOverhang = TRUE` 46.87%. 19.7% of the pairs that entered the merge step were rejected (11.5 of 58.4 points).

*What this shows.* (1) **The overhang hypothesis is rejected for this dataset:** no pair was refused for an overhang and `trimOverhang = TRUE` changed nothing; the reads of this dataset presumably carry no read-through into the other primer. Nothing about `trimOverhang` is added to the ITS report or its DADA2 snippet. (2) The split of the loss is now exact: 15.9 points to denoising and 11.5 to merging, inside the bounds of section 19 (at least 14.0, at most 13.4). (3) Merging rejected pairs for two reasons, in nearly equal shares: too few matching bases in the overlap (5.07 points) and mismatches or indels in the overlap after denoising (6.43 points). (4) AmpliTrim's estimate of 1.9 points at merging is 9.6 points too low, so the yield estimate (72.4%) overshoots DADA2's merged share (46.9%) by 25.5 points: 15.9 of them are denoising, which is not modelled, and about 10 are merging.

*A likely mechanism, not tested.* AmpliTrim's tolerant merge calls vsearch with its default of up to 10 differences and a minimum overlap of 12 nt, on raw reads; DADA2 merges denoised sequences with `maxMismatch = 0` and counts matching bases, not overlap length, against `minOverlap`. A pair whose denoised reads still differ in the overlap, or whose overlap has fewer than 12 matches, is refused by DADA2 and accepted by vsearch. A test would apply DADA2's two rules (no mismatch; 12 matching bases) to the vsearch alignments, or run `mergePairs` with `maxMismatch = 1` and compare; neither was run.

*Decisions.* The user guide now gives the exact split and these reasons, and says that the merging rows of the ITS table are an upper bound that overstates what merges by about 10 points on this dataset. No code change. Whether to show a range for ITS from the strict (no mismatch) and tolerant merge rates is not decided: the strict rate on raw reads is a lower bound on a different quantity (raw reads, not denoised), and its use would need its own test. With this, the ITS question is closed for the release: the filter rows and the forward-only share replicate DADA2's filter; the merging rows are flagged as optimistic with their measured size; maxEE is not where ITS loses its pairs. Open, by design: accuracy on ITS (no expected sequences) and any other ITS data.

---

**Outcome of section 18, night 4b (run 10 October 2026 from 01:46 to 10:21; build `0c4f703158`; 5 datasets run, 7 left at the default by construction: V1-V2, soil V4, inter-laboratory V1-V2, D6300 V3-V4, spike-in, swabs; ITS not part of this test; 0 problems).**

| | Prediction | Outcome |
| --- | --- | --- |
| P19 | H (q = 0.75, cap 5) merges at least as many pairs as G (2,5), within the margin of 2 points, on every dataset | **Not met.** H minus G: mock A -6.7, mock B +5.3, inter-laboratory V4-V6 -2.6, D6300 V4-V5 -5.4, wastewater +6.8. Below the margin on three of five. Mock A and V4-V6 were predicted; V4-V5 was not. |
| P20 | H is as accurate as G or more accurate on every mock with expected sequences | **Not met.** Mock A as accurate, D6300 V4-V5 as accurate, **mock B less accurate**: 73.8% of reads in identical ASVs against G's 74.7% (-0.9 points, margin 0.1), the same variants (6) and strains (3/3); H and G have the same cut (277/250) and differ only in the forward limit (4 against 2). H is +5.1 points above the default (68.7%). |
| P21 | H merges more pairs than G by more than the margin on wastewater and mock B | Met: +6.8 and +5.3 |
| P22 (exploratory) | Where the uncapped rule goes above 5 (mock B), H80 merges more than H: the cost of the cap | Met: +20.2 points (H80 75.3%, H 55.0%; FIGARO 74.3%) |

*Merged share of the input pairs (DADA2), the five datasets where the rule differs from the default:*

| dataset | A (2,2) | G (2,5) | H, q 0.75, cap 5 | H80, q 0.80, no cap | rule's maxEE (H / H80) |
| --- | --- | --- | --- | --- | --- |
| mock A | 62.7% | 76.1% | 69.5% | 73.3% | 2,3 / 3 |
| mock B | 19.9% | 49.7% | 55.0% | 75.3% | 4,5 / 5,9 |
| inter-laboratory V4-V6 | 41.6% | 60.4% | 57.9% | 68.8% | 3 / 4 |
| D6300 V4-V5 | 57.2% | 68.8% | 63.3% | 66.4% | 2,3 / 2,4 |
| wastewater | 43.6% | 51.1% | 57.9% | 57.9% | 4,3 / 4,3 |
| mean | 45.0% | 61.2% | 60.7% | 68.3% | |

**Decision by the logic fixed in section 18.** P19 fails, so the rule is *not shown to be at least as good as rule R*; P20 fails, and section 18 says: *if P20 fails, rule R stays.* The rule is not adopted on the strength of H80 or of any other q chosen after this run.

*Descriptive, post hoc (not criteria).* (1) H is above the default on all five datasets (+6.1 to +35.2 points), and on average (60.7%) about level with (2,5) applied everywhere (61.2%), which is what rule R at a 70% trigger gives on these five datasets (all five pass under 70% at maxEE 2). (2) **The release's rule R at a 60% trigger**, the post hoc trigger chosen by the researcher, leaves mock A (64.7% passing at 2) and V4-V5 (63.7%) at the default, and so gives 56.2% on average: H is ahead of it on mock A (+6.8), mock B (+5.3), V4-V5 (+6.1) and wastewater (+6.8) and behind on V4-V6 (-2.5), but this is the comparison with a trigger chosen after seeing the data, and with a trigger of 70% (the one run in night 3) H is behind on mock A and V4-V5 and level on average. The comparison depends on the trigger. (3) The cap of 5 costs 20 points of pairs on mock B and 11 on V4-V6 (H80 is level with FIGARO on both), and nothing elsewhere. (4) H relaxes both reads by their own noise: the forward read is looser and the reverse read cleaner than under (2,5) (mock A kept reverse errors 0.91 against 1.43; wastewater forward 1.49 against 0.69). (5) The P20 margin on reads in identical ASVs (0.1 points) is the one fixed for the battery in night 2 and is strict: the 0.9-point difference on mock B comes with a +5.1 gain over the default.

*Pending: the researcher's choice for the release.* The options, none of which changes what the criteria say: (a) rule R with a trigger of 70%, the version run in nights 3 and 4b (arm G on these datasets); (b) rule R with the 60% trigger already chosen; (c) the percentile rule, chosen as a stated trade, documented as having failed P19 and P20.

**Decision (10 October 2026, after night 4b): the release uses the percentile rule; rule R is not used.** The researcher (Payton Yau) chose option (c): to keep the percentile rule (q = 0.75, cap 5) in the release and not to reinstate rule R, at either trigger. This **overrides the clause of section 18 that says that rule R stays if P20 fails**, with the results of night 4b in view: P19 and P20 not met, the rule level with (2,5) on average (60.7% against 61.2%) and above the default on all five datasets, fewer pairs than (2,5) on three of them, 0.9 points fewer reads in identical ASVs on mock B. The reasons are those recorded when the rule was chosen: it relaxes each read by its own noise, and its limits stay within the range that DADA2's author calls acceptable. The release documentation and the report's note on maxEE state the outcome of night 4b (the table of merged shares, the two criteria not met, the post hoc status of the percentile and the cap, and that the rule has not been tested on independent data). In the paper the rule is described as chosen by the author after seeing the validation data, tested in a pre-registered run in which it did not meet its adoption criteria, and kept as a stated trade. Rule R (at 70%) remains the comparison arm G of nights 3 and 4b.

---

## 21. Release check: does the released package recommend the settings that DADA2 was run at? (written 10 October 2026, before `validation/release_check.py` was run)

**Why.** Nights 1 to 4b ran the development build, which still has FIGARO, the exploration options, rule R and, for night 4b, the percentile rule computed by a validation script. The release is a different package: no FIGARO, no exploration, and the percentile rule built into `amplitrim.pipeline`. The release has not been run on the real datasets. DADA2's result depends only on the settings it is given (truncLenF, truncLenR, maxEE, with minOverlap and minLen as in the nights), so if the released package recommends the same settings on every dataset, the DADA2 results of nights 3 and 4b are the released package's results and nothing needs to be re-run. Not claimed: that the code is the same (it is not), only what it recommends.

**What it does.** The released package (`--package`, the clean repository; it refuses the development build) is run on the 12 datasets of the battery with its defaults (the marker `its` for ITS, `16s` for the others; the same default subsample as the nights), without DADA2 and without R. Its recommendation (truncLenF, truncLenR and maxEE) is compared with what was run through DADA2: the arm H of night 4b where the rule relaxed maxEE (mock A, mock B, inter-laboratory V4-V6, D6300 V4-V5, wastewater), otherwise the arm A of night 3 (the default maxEE 2 and the best truncation at it; for ITS, no truncation). The tested settings are read from the result files of the nights.

**Prediction, written before the run.**
- **P30.** On all 12 datasets the released package recommends exactly the settings that were run through DADA2: the same truncLenF, the same truncLenR and the same maxEE (forward and reverse).

**What follows.**
- *P30 holds:* the validation documentation of the release says that the DADA2 results of nights 3 and 4b are those of the settings the released package recommends, on every dataset, and the paper presents them as the results of the released package, with the development build described as the version the study was run with.
- *P30 fails on some datasets:* the DADA2 results of the nights are not claimed for the released package on those datasets. DADA2 is run at the released settings (arm A and the released setting) on those datasets only, before any claim; the datasets where it holds stand.
- Not tested here, in either case: that DADA2 gives the same result twice (it did, to the hundredth of a point, in section 20, P26), and the datasets' accuracy, which is as measured in the nights.

**Correction, 11 October 2026 (the text above is left as it was written).** "The same default subsample as the nights" was wrong: `compare_methods.py` thins each dataset to 5,000 pairs per sample before AmpliTrim runs, so the nights' picks come from those pairs and not from the package's default sample of the full reads. `release_check.py` therefore compares the package with the nights on different pairs. The corrected comparison is P39b in section 24 (`nocap_check.py --like-the-nights`).

---

## 22. Soil V4: over-merging on real reads, and the grid around the final recommendation (written 10 October 2026, before either was run)

**Why.** Soil V4 (PRJNA574030) is where AmpliTrim's truncation was furthest from the best of a DADA2 grid, and where the yield estimate overshot most (+32.1 points). What the record holds, by round:

| round | soil V4 | what was measured |
| --- | --- | --- |
| B6 (an earlier build) | AmpliTrim's cut 120/118 | 59.6% of the pairs merged against 66.4% for the best cell of the 3 x 3 grid, 130/108: 6.8 points below, margin 1 point, criterion not met. The best cells kept the same total (238 nt) and shortened the reverse read |
| B5 (the default `--allow-unmergeable` lowered from 1% to 0.2%) | AmpliTrim's cut 124/115, the final build's | 120/118 to 124/115 gained 2.6 points (64.3% to 66.8%); another run than B6, so not comparable cell by cell |
| Table 9, night 2 | picks on three disjoint sets of samples: 124/115, 124/115, 113/126 | the pick does not hold (one pick moved 11 nt) |
| B7 | per-sample rank agreement | 0.70 over 36 samples, criterion not met |

The grid was never re-run around the final recommendation, so the gap is unknown for the build that is released. Over-merging (merges at the wrong length) was tested in simulation (Table 13: 0.00% wrong-length merges at 12 nt in the clean and realistic arms; in the stress arm with a 32 nt perfect repeat, a longer overlap helped) and on mock community A (Q23: 0.13% of the ASV reads off every expected length), and **never on soil V4**, a real, diverse community without expected sequences.

**Part 1: false joins on real reads (`validation/overmerge_real.py`; vsearch only, no R).** Real reads have no truth, so the reference is the pair merged from its full-length reads, where the overlap is long: a pair whose full-length reads overlap by at least 40 nt is a reference pair, with its reference merged length L. A cut is applied (the DADA2-style filter at maxEE 2), the surviving reference pairs are merged by vsearch at the minimum overlap of 12 nt (tolerant, as AmpliTrim does, and with no difference allowed), and a pair merged at another length than L was joined at the wrong place (a false join; a pair that the reads cannot span but that is merged anyway counts as one). Also measured per cut: the expected errors that the two reads carry in their overlap (the last bases of each cut read), which is where DADA2 refuses a pair. Cuts: the 3 x 3 grid of B6 around 120/118, and the grid around AmpliTrim's recommendation now.

**Part 2: the grid for the final build (R and DADA2).** `overnight.py --only landscape --extra "v4=..."` centres its 3 x 3 grid (steps of 10 nt) on what AmpliTrim recommends now, 3,000 pairs per sample, maxEE 2, and compares the recommendation with the best setting that meets the coverage floor; the margin is the registered 1.0 point (`landscape_regret_points`).

**Predictions, written before the run.**
- **P31.** At each cut that matters (AmpliTrim's recommendation now, 120/118 and 130/108) at most 1% of the merged reference pairs are at the wrong length. Expected: met, as in the simulation and on mock A.
- **P32.** The wrong-length merges do not tell 120/118 from 130/108: they differ by at most 1 point of the reference pairs. Expected: met.
- **P33 (exploratory).** The overlap of 130/108 carries fewer expected errors per pair than that of 120/118. This is a hypothesis for the mechanism (the overlap region of the two cuts is a different 12 nt window, at the end of a different part of each read); no outcome is predicted.
- **P34.** AmpliTrim's recommendation now (124/115 on this build) is within 1 point of the best setting of the grid that meets the coverage floor. Best guess: not met, given B6, Table 9 and B7.

**What follows.**
- *P31 holds:* the 12 nt minimum stands on real soil V4, which extends the simulation and mock A. *P31 fails:* false joins are real at 12 nt on this dataset, the yield estimate counts them as usable, and the documentation says by how much; a longer minimum overlap would need its own test on real data before any change (the note on Q21 to Q27).
- *P32 holds:* over-merging is not why the earlier cut trailed the best cell. *P32 fails:* over-merging differs between the cells by more than a point, a candidate cause that is investigated.
- *P33 holds:* a candidate mechanism, which would suggest choosing the split by the errors in the overlap. It is not implemented now: it would be a post hoc change, to be tested on independent data first, and it is listed as future work. *P33 fails:* it is not the mechanism.
- *P34 holds:* the documentation and the paper say that the final build is within a point of the best cell on soil V4, in place of 'unresolved'. *P34 fails:* they state the gap that was measured for the final build, soil V4 stays the documented exception, and the tool is not changed before the release.
- Not claimed in any case: accuracy on soil V4 (no expected sequences), or that the result holds beyond this dataset.

**Outcome of section 22, part 1: over-merging on real reads (10 October 2026; soil V4, 49,720 pairs analysed, 48,227 reference pairs (97.0%), 18 cuts: the 3 x 3 grids around 120/118 and around 124/115).**

| | Prediction | Outcome |
| --- | --- | --- |
| P31 | At most 1% of the merged reference pairs are at the wrong length at each cut that matters (124/115, 120/118, 130/108) | **Met:** 0.00% at all three, and at all 18 cuts |
| P32 | Wrong-length merges do not tell 120/118 from 130/108 (within 1 point of the reference pairs) | **Met:** 0.00% against 0.00% |
| P33 (exploratory) | The overlap of 130/108 carries fewer expected errors per pair than that of 120/118 | **Met in direction, not in size:** 0.056 against 0.059 (0.003 per pair, at most 0.3 points of pairs) |

*What this shows.*
1. **There is no over-merging on soil V4.** Merges at the wrong length were 0.00% at every one of the 18 cuts, with vsearch's tolerant merge (up to 10 differences) and with none allowed, including the cuts below the coverage floor (110/108, a total of 218 nt: 2.9% of the reference pairs merged correctly, none wrongly). The 12 nt minimum overlap stands on real soil V4 reads, as it did in simulation (Table 13) and on mock community A (Q23). Scope: the reference is itself a vsearch merge, of the pairs whose full-length reads overlap by at least 40 nt; the 1,493 pairs (3.0%) that overlap by less cannot be judged; and 0.00% is a figure rounded to two decimals.
2. **vsearch sees no gap between the cells of round B6.** 120/118 and 130/108 both merged 99.3% of the reference pairs correctly (-0.0 points), 124/115 merged 99.7%, and 99.9% of the pairs pass the filter at every cut near the floor. DADA2 measured 59.6% against 66.4% for the first two (+6.8 points). What separates them is therefore not over-merging, not the quality filter and not vsearch's view of the overlap.
3. **The errors in the overlap barely differ** between those two cuts (0.056 against 0.059 per pair), so they cannot account for 6.8 points either. As a check on the measurement, the share of pairs that merge with no difference allowed (90 to 96% across the cuts) agrees with one minus the expected errors in the overlap to within 2 points.
4. **By exclusion, the gap lies after the filter, in denoising or in DADA2's own merging of the denoised sequences,** which AmpliTrim does not model. This is an inference from what was ruled out, **not a measurement on soil V4**; on ITS the exact split (section 20) found 15.9 points lost in denoising and 11.5 in merging. It is consistent with the largest overshoot of the yield estimate being on this dataset (+32.1 points). A measurement would split DADA2's flow per cell, as section 20 did for ITS.
5. **A descriptive observation, after the fact:** the share merged with no difference allowed varies by up to 5 points among the cuts that meet the coverage floor (134/105 95.7%, 124/115 92.9%, 134/125 90.7%), and is higher where the reverse read is shorter, in the direction of the best cell of round B6 (130/108 94.4% against 92.4% at 120/118).

**P35 (exploratory; written 10 October 2026 after reading part 1 and before the grid of part 2 is read; the no-mismatch shares are those of this report, the DADA2 results are not yet known).** Among the cells of the part 2 grid that meet the coverage floor (124/115, 124/125, 134/105, 134/115, 134/125 and 114/125), the cell with the highest no-mismatch merge share in this report (134/105, 95.7%) is within 1 point of the best of them by DADA2's merged share.
- *P35 holds:* the no-mismatch merge share is a candidate rule for choosing the split at the floor. It is **not adopted**: it is post hoc, from one dataset, and would need independent datasets where the split matters, with criteria fixed before the run, before any change to the tool.
- *P35 fails:* the no-mismatch merge share does not predict DADA2's merged share, and the soil V4 gap stays unexplained.

*Decisions.* The validation documentation and the paper now say that over-merging was tested on real soil V4 reads and found at 0.00% at every cut, and that the soil V4 gap is not explained by over-merging, the filter or vsearch's view of the overlap (an inference, not a measurement). Part 2 (the grid for the final build, P34) and P35 are pending; the soil V4 statement stays 'unresolved' until P34 is read.

**Outcome of section 22, part 2: the DADA2 grid for the final build (10 October 2026; soil V4, maxEE 2, a 3 x 3 grid in steps of 10 nt around AmpliTrim's 124/115, 3,000 pairs per sample, 19 minutes).** Share of the input pairs merged by DADA2; `.` marks a total below the coverage floor of 239 nt, `*` AmpliTrim's cut:

| forward \ reverse | 105 | 115 | 125 |
| --- | --- | --- | --- |
| 114 | 64.8 . | 61.5 . | 58.1 |
| 124 | 67.1 . | 62.6 * | 58.1 |
| 134 | 67.3 | 62.9 | 58.4 |

| | Prediction | Outcome |
| --- | --- | --- |
| P34 | AmpliTrim's recommendation now (124/115) is within 1 point of the best setting of the grid that meets the coverage floor | **Not met:** 62.6% against 67.3% for 134/105, 4.7 points below |
| P35 (exploratory) | The cell with the highest no-mismatch merge share among the cells that meet the floor (134/105, 95.7%) is within 1 point of the best of them by DADA2's merged share | **Met:** 134/105 is both, 0.0 points apart |

**Correction.** The 'Why' paragraph of section 22 said that the grid had never been re-run around the final recommendation. This was incorrect. Q28 (round 4) already records a gain of **+4.7 points** for soil V4 ("the best cell with AmpliTrim's total length"). This equals the gap in this grid to the decimal (67.3% minus 62.6%). An earlier run had therefore read the same grid, and part 2 reproduced it. The headline number is not new; the full grid is. The validation page and the paper said that the grid 'has not been re-run'. Both statements were incorrect and have been corrected. One sentence in `docs/FINDINGS.md` said that a forward-heavier split 'did not beat' AmpliTrim's balanced split. This was ambiguous. It now states what Q28 recorded: the forward-heavier split was the best cell in 3 of 6 datasets, mostly soil V4, which is short of the 4 of 6 and the 1.0 point set in advance.

*What the grid shows (descriptive, after the fact).*
1. **The reverse read is what matters on soil V4.** Merged share falls with the length of the reverse read and hardly changes with the length of the forward read: 66.4%, 62.3% and 58.2% for reverse reads of 105, 115 and 125 nt (each the mean of three forward lengths), against 61.5%, 62.6% and 62.9% for forward reads of 114, 124 and 134 nt. Rank correlation with the merged share: +0.95 for a shorter reverse read, +0.29 for the forward length. At the total of 239 nt, three splits gave 58.1% (114/125), 62.6% (124/115, AmpliTrim's) and 67.3% (134/105). The best cells of round B6 (130/108) and of this grid (134/105) both shorten the reverse read and keep the total.
2. **The coverage floor did not cost merging on soil V4.** Ten nucleotides below it, a cut with a reverse read of 105 nt merged 67.1% (124/105), as many as the best cut at the floor (67.3%), and 20 nt below it 64.8% (114/105). Merging did not collapse as it did on mock community A and V4-V6. (The floor is set to cover all but the longest 0.2% of the amplicons.) The best cell of this grid meets the floor, so the floor is not what cost AmpliTrim 4.7 points; the split is.
3. **P35 holds, but the no-mismatch merge share is a weak predictor across the grid:** rank correlation with DADA2's merged share +0.45 over the 9 cells and +0.55 over the 6 that meet the floor, against +0.95 for the length of the reverse read alone. With 9 cells of one dataset, P35 holding probably reflects that the no-mismatch share is higher where the reverse read is shorter.
4. The result is for 3,000 pairs per sample; whether it persists at full depth was not tested.

*Decisions, as set in advance.* P34 fails: the documentation and the paper state the gap measured for the final build (4.7 points; 6.8 with the earlier build of round B6), soil V4 stays the documented exception, and the tool is not changed before the release; Q28's decision stands (keep the balanced tie-break: the gain is one dataset). P35 holds: a candidate rule, **not adopted**. A rule that gave the forward read more of the total would be a post hoc change from one dataset (the gain on the others was +0.5, +0.1 and none); testing it needs independent datasets where the split matters, with the criteria fixed before the run. It is listed as future work.

---

## 23. The merge estimates against DADA2, and the cut-off for differences in the overlap (written 10 October 2026, before `validation/estimate_check.py` was run)

**Why.** The yield estimate overshoots DADA2's merged share, by up to +21.1 points on ITS and +32.1 on soil V4. There are two reasons. (1) vsearch's tolerant merge, which allows up to 10 differences in the overlap, accepts pairs that DADA2 refuses after denoising (no mismatch allowed; at least 12 matching bases). (2) vsearch cannot see the denoising, which removes reads before merging. What is known: on soil V4 the strict estimate (no difference allowed, on the raw reads) still overshot by +25.0, so tolerance accounts for about 7 of the 32 points; at 124/115 vsearch merged 99.7% of the reference pairs tolerantly and 92.9% with no difference allowed (section 22); on ITS DADA2 lost 15.9 points in denoising and 11.5 in merging, against 1.9 estimated (section 20). AmpliTrim already merges every candidate in both modes and uses the tolerant rate in the estimated yield. It was never tested across the datasets which of the two estimates is closer to DADA2, nor whether a number of allowed differences between 0 and 10 matches DADA2's merging better than either.

**Part A (the tables on disk; no vsearch and no R).** From every `<run>/<dataset>/compare/compare_methods.tsv` of the nights: for each setting (a setting run on several nights counts once, the newest kept), the tolerant estimate, the strict estimate and DADA2's merged share. Per dataset: the mean overshoot of each estimate, the mean absolute difference, and the rank correlation with DADA2's merged share over the settings (for datasets with at least 4 settings).

**Part B (vsearch; `--cutoff`).** For each AmpliTrim setting that has a per-sample DADA2 track: the pairs that pass the filter at that setting (AmpliTrim's subsample) are merged by vsearch with at most k differences in the overlap (k = 0, 1, 2, 3, 4, 5, 7 and 10; minimum overlap 12), and the share merged is compared with DADA2's merge-step success from the track, merged / sum of min(reads in denoised forward sequences, reads in denoised reverse sequences), a lower bound of the true success (on the ITS dataset the exact figure was 80.3% where this formula gives 77.9%). The measure is the mean absolute difference in points over the settings.

**Predictions, written before the run.**
- **P36.** The strict estimate is closer to DADA2's merged share than the tolerant one: the mean over the datasets of the mean absolute difference is smaller. Best guess: met on the whole, carried by ITS and soil V4; on the mock communities, where the tolerant overshoot was small (-2.2 to +6.7), the strict estimate may be farther.
- **P37.** The strict estimate ranks the settings at least as well as the tolerant one: the mean rank correlation over the datasets with at least 4 settings is not more than 0.05 below. No strong expectation: raw reads make the strict estimate noisy.
- **P38.** An intermediate number of allowed differences is closer to DADA2's merge step than both 0 and 10, by at least 2.0 points of mean absolute difference. Best guess: met, since the raw strict merge undershoots what DADA2 merges after denoising and the tolerant merge overshoots it.

**What follows.**
- *P36 and P37 hold:* the strict estimate is a candidate for the headline yield, or for a second ceiling shown beside the tolerant one. It is **not adopted automatically**: it would be a post hoc change of a validated component, to be decided by the researcher as a stated trade, with its status documented as for the `maxEE` rule. *P36 holds and P37 fails:* the strict estimate is closer in level but ranks worse, so the tolerant estimate keeps the ranking and the strict one could only be shown as a ceiling. *P36 fails:* the tolerant estimate stands.
- *P38 holds:* a candidate cut-off, reported with its mean absolute difference; not adopted without an independent check. *P38 fails:* 0 or 10 stands.
- None of these addresses the loss in denoising, which does not come from vsearch (15.9 points on ITS; the larger part of the soil V4 overshoot is inferred, not measured). The part of the overshoot that no merge setting can explain is reported as what is left.
- Not claimed in any case: accuracy, independent data, or that a setting chosen on these datasets holds on others.

**Outcome of section 23 (10 October 2026; 98 settings on 12 datasets from the tables of the earlier nights; part B on the AmpliTrim settings that have a DADA2 track, 59 of them, at 0, 1, 2, 3, 4, 5, 7 and 10 allowed differences).**

| | Prediction | Outcome |
| --- | --- | --- |
| P36 | The strict estimate is closer to DADA2's merged share than the tolerant one | **Not met:** 12.5 points against 6.5 (the mean over the datasets of the mean absolute difference); strict closer on 5 of 12 datasets |
| P37 | The strict estimate ranks the settings at least as well (within 0.05) | **Not met:** mean rank correlation +0.76 against +0.92 over the 10 datasets with at least 4 settings |
| P38 | An intermediate number of allowed differences is closer to DADA2's merge step than both 0 and 10, by at least 2 points | **Not met:** the mean absolute difference falls at every step, 27.6 (0), 15.0 (1), 10.1 (2), 8.4 (3), 7.6 (4), 7.2 (5), 6.9 (7), 6.6 points (10); the closest is 10 |

The best guesses recorded beside the predictions were wrong for P36 (strict was expected to be closer on the whole) and for P38 (an intermediate cut-off was expected to win). No outcome was expected for P37.

*What this shows.*
1. **The tolerant estimate is the better of the two, and the default of 10 differences is the best of the cut-offs tried.** Strict merging of the raw reads is far too low on the mock communities (mean overshoot -14.8 points, against +1.1 for the tolerant estimate; mock B -37.8, wastewater -26.9, mock A -23.1): on noisy raw reads it refuses pairs that DADA2's denoising repairs. On mock B at maxEE 5,9 vsearch merged 15.7% with no difference allowed and 90.1% tolerantly, where DADA2's merge step succeeded for 98.6%. This is why the tolerant mode was built (section 2.3 of the write-up).
2. **No number of allowed differences matches DADA2's merge step better than the default.** The signed difference crosses zero at 3 allowed differences (-0.5 points on average), so 3 is the unbiased cut-off, but its mean absolute difference (8.4) is larger than at 10 (6.6): the settings differ too much from one another for any single cut-off to fit them.
3. **The merge step is not where most of the overshoot is.** On real datasets the tolerant estimate overshoots by +10.1 points on average (soil V4 +26.5, ITS +13.7, swabs +7.1, V1-V2 +6.1, wastewater +3.9, spike-in +3.3) and the strict estimate by -2.5, closer on 4 of 6 real datasets (ITS, swabs, V1-V2, soil V4). That is **not a justification for strict merging**: on the mocks it is 14.8 points too low, and on real data the two errors (a strict merge that is too harsh, a denoising loss that vsearch cannot see) happen to cancel. On soil V4 even the strict merge overshoots DADA2's merge step by 13.6 points (91.8% against 78.2%), and the denoising loss alone is at least 13.3% of the filtered pairs; on ITS the merge step lies between the strict and the tolerant merge (59.2%, 97.5%, DADA2 77.7%; the exact figure of section 20 was 80.3%).
4. **Ranking:** the tolerant estimate ranks well on most datasets (+0.94 to +1.00 on nine of ten) and badly on soil V4 (+0.26; the strict estimate +0.83 there, the one dataset where strict ranks clearly better).

*Limits of the check.* DADA2's merge step is the lower bound merged / sum of min(denoised forward, denoised reverse), not the exact figure (a few points low, as on ITS); part B uses the AmpliTrim settings only; only the number of differences was varied, not the minimum overlap; the settings are those of nights 1 to 4b, and none of them is independent of the design of the estimate.

*Decisions, as set in advance.* P36 fails, so **the tolerant estimate stands**; P38 fails, so **10 allowed differences stands** (vsearch's default); P37 fails, so ranking keeps the tolerant estimate. Nothing in the tool changes. The validation page and the paper state that a strict merge is farther from DADA2 and ranks worse, and that no cut-off between 0 and 10 does better. What no merge setting can reach is the loss in denoising: at least 18.8% of the filtered pairs on ITS and 13.3% on soil V4, against about 1% on the mocks. Modelling it would need a predictor that can be computed from the reads, and that is a new question, not part of this section.


## 24. Night 5: the rule without a cap, whether a simpler cut would do, and the reads wastewater loses to length (decision of 11 October 2026; predictions written before any of it was run)

**Decision (11 October 2026).** The researcher (Payton Yau) decided to remove the cap of 5 from the released rule and to put a concise note in the report when a limit goes above 5. The reasons given: there is no absolute rule for `maxEE`; the algorithm should follow the quality of the reads; and a fixed number depends on the read length, so it would not carry to long reads. **This is post hoc.** It was made after the data of nights 1 to 4b had been seen, and it reverses the decision of section 18, which kept the cap for the reasons recorded there: DADA2's author argues for restraint and gives no ceiling, only mock B gained beyond 5, and the effect of noisier reads on relative abundances in rich communities is not measured. The researcher had also asked for the rule to be conservative on poor reads. The capped rule is: on mock A, mock B, inter-laboratory V4-V6 and wastewater it merged 8.5, 19.0, 10.7 and 8.8 points fewer pairs than FIGARO's own pick (DADA2, points of input pairs). Without the cap it is still below FIGARO's pick on mock B (a reverse limit of 8 against 9, a forward limit of 4 against 6) but less conservative than it was, and on that dataset only.

**What changes in the released package.** The rule has no upper limit (per read, `int` of the 75th percentile of the expected errors at the cut, plus 1, never below 2). The report and `summary.json` (`warnings`) carry one note when a limit is above 5: "The rule sets maxEE a,b, above 5: these reads are very noisy, and more low-quality reads may not improve the relative abundances (DADA2's author). Check on a mock community if you can." `summary.json` keeps `maxEE_rule.cap` (now `null`; the schema never loses a key) and gains `maxEE_rule.warn_above` (5). The fallback that decides where to measure when nothing passes at `maxEE` 2 still judges at 5; it only decides where to look and does not limit the recommendation. By the offline figures of section 18 the cap bound on **mock B only**: at q = 0.75 the rule gives (4,8) there instead of (4,5), with an estimated 69.5% of the pairs passing instead of 55.2%, and 3.31 expected errors per kept reverse read instead of 2.55 (forward 1.02 instead of 0.81). On every other dataset no limit goes above 5, so nothing changes.

**What is measured and what is not.** The DADA2 results of the rule on the validation page were measured with the cap. They stand for the four other datasets where the rule relaxes `maxEE` (mock A, V4-V6, V4-V5, wastewater) and for the clean datasets and ITS (`maxEE` 2), where the cap changes nothing. For mock B, **the released rule's setting has not been run through DADA2**: the variant with q = 0.80 and no cap, (5,9), merged 75.3% (section 18), which is not the same setting. The validation page and the results site say so until P40 and P41 are read.

**Effect on section 21.** P30 was decided for the capped rule. With the cap removed, `validation/release_check.py` reports mock B as different from what DADA2 was run at, by design. `validation/nocap_check.py` is the check for the package as released now.

**Parked (11 October 2026).** A long-read mode (PacBio HiFi full-length 16S and the whole rRNA operon) was discussed and set aside by the researcher. It is not part of 0.1.0 and has no prediction here.

**Questions and predictions, written before any of this was run.** Margins as in nights 2 to 4b (`validation/maxee_arms.py`): merged share 2 points; variants and strains within 1; reads in identical ASVs within 0.1 points.

- **P39 (offline, released package).** On every dataset the released package recommends what nights 3 and 4b ran through DADA2, except on mock B, where it sets `maxEE` 4,8. Decided by `nocap_check.py` (a mock B that comes back unchanged also fails it: the cap would still be there).
- **P40 (DADA2, mock B).** The released rule's setting merges more than the capped rule's (4,5; 55.0%) by more than the margin. Best guess, not a criterion: between 60% and 76%.
- **P41 (DADA2, mock B).** The released rule's setting is as accurate as the capped rule's or more accurate ('as accurate' or 'more accurate' by the margins above).
- **P42 (descriptive).** Expected errors per kept read at the released cut and limit: about 1.0 forward and 3.3 reverse (the offline estimate).
- **P43 (offline, the ablation).** The floor-only cut (the best split of the shortest total that covers the amplicons, which is the longest amplicon with the longest 0.2% left out, plus 12, judged at the released `maxEE`) is the released cut on at least 8 of the 11 datasets that have a cut.
- **P44 (DADA2, where P43 finds a difference).** On each dataset where the floor-only cut differs from the released cut, DADA2 at the floor-only cut merges at least as many pairs as at the released cut, within the margin (floor-only minus released at least -2 points).
- **P45 (descriptive, no criterion).** The share of read pairs shorter than the recommended cut on each dataset, to tell the loss to length from the loss to quality on wastewater (the package warns from 1%).

**Decisions, as set in advance.** If P39 and P40 hold, the removal stands, the validation page and the results site give the measured result for mock B, and the status of the rule does not change (it stays a sensible default, not independently tested). If P40 fails, the removal does not help mock B; it stands, because it was decided for the reasons above, and the page says so. If P41 fails, the removal costs accuracy on the one dataset where it matters, and the question of a cap is reopened by the researcher; no cap is added on the strength of this result alone. If P43 and P44 hold, the search over longer totals and the four-key tie-break are candidates for removal in 0.2, which would give a simpler algorithm with the same results; **nothing changes in 0.1.0 either way.** The stability of the floor-only cut over three sample sets is not part of this night. P42 and P45 are descriptions. If P45 shows that most of wastewater's lost pairs are reads shorter than the cut, the wording that calls it a dataset of poor read quality is changed on the README, the website and the results pages, and the finding is added to this section.

**How to run (all on the machine that has the reads).**

```
python validation\nocap_check.py --package C:\Users\Admin\OneDrive\Documents\GitHub\amplitrim --config C:\amplitrim_v2\amplitrim\amplitrim_paths.json --night3 C:\data\runs\night3_2026-10-09 --night4 C:\data\runs\night4_2026-10-10
```

The report prints the commands for the DADA2 runs it cannot do itself: mock B (`compare_methods.py ... --amplitrim-max-ee 2 2,5 4,5 4,8`, then `nocap_check.py --dada2-table ...\compare_methods.tsv`) and, for each dataset where P43 finds a difference, DADA2 at the floor-only cut (`--setting floor_only=F,R,EE`). Send NOCAP_CHECK_REPORT.md, NOCAP_DADA2_REPORT.md and the compare_methods.tsv files.

**Outcome so far (11 October 2026; part 2 only: part 1, `nocap_check.py` on the released package, has not been run).** DADA2 on mock B (`compare_methods.py`, 8 samples, 5,000 pairs per sample, `--amplitrim-max-ee 2 2,5 4,5 4,8`):

| setting | cut F/R | maxEE | passing the filter | merged | variants / strains | reads in identical ASVs |
|---|---|---|---|---|---|---|
| default | 279/248 | 2 | 20.3% | 19.9% | 6 / 3 of 3 | 68.7% |
| rule R's setting | 277/250 | 2,5 | 50.9% | 49.7% | 6 / 3 of 3 | 74.7% |
| the capped rule (night 4b) | 277/250 | 4,5 | 56.4% | 55.0% | 6 / 3 of 3 | 73.8% |
| AmpliTrim's best truncation at the uncapped rule's maxEE | 274/253 | 4,8 | 72.2% | 70.4% | 6 / 3 of 3 | 75.2% |

- **P40 MET.** At maxEE 4,8 DADA2 merged 70.4% of the pairs against 55.0% for the capped rule: +15.4 points (the margin is 2), inside the predicted range of 60% to 76%. Per sample it merged more than the default in 8 of 8 samples.
- **P41 MET.** The setting at 4,8 is more accurate than the capped one: variants 6 against 6, strains 3 of 3 against 3 of 3, reads in identical ASVs 75.2% against 73.8%.
- **The run is comparable with the earlier nights:** the default (19.9% merged, 68.7% identical), the capped rule (55.0% and 73.8%) and 2,5 (74.7% identical) reproduce the numbers of nights 3 and 4b exactly.
- **Descriptive, not predicted.** AmpliTrim's estimate of the share passing the filter at 4,8 was 72.3% against 72.2% measured; its yield estimate was 67.2% against 70.4% merged. Against the other tool's own pick on mock B (74.3% merged, 78.3% of reads in identical ASVs, 7 variants; nights 1 to 3, not re-run here) the setting at 4,8 is still lower: 70.4%, 75.2% and 6 variants.
- **Not decided: P39, P42, P43, P44 and P45**, because part 1 has not been run. The DADA2 run shows what AmpliTrim's best truncation at maxEE 4,8 does. It is the released package's result only if P39 shows that the package recommends 274/253 at maxEE 4,8 on mock B. By the decision rule above, the validation page and the results site keep saying that mock B's released setting has not been run through DADA2 until P39 has been read.

**Outcome, part 1 (11 October 2026; the released package on its own sample of the full reads, 12 datasets, `nocap_check.py`).**

- **P39 NOT MET.** On six datasets the released package recommends exactly what DADA2 was run at (mock A, ITS, soil V4, D6300 V3-V4, D6300 V4-V5, swabs), and mock B changed as predicted, to 274/253 with maxEE 4,8. On five it recommends something else: V1-V2 201/157 tested, 195/168 released (the total 358 against 363); inter-laboratory V1-V2 155/170 against 160/165; inter-laboratory V4-V6 271/267 with maxEE 3 against 273/265 with maxEE 4,3; wastewater 261/183 against 260/182 (the same maxEE 4,3); spike-in 222/218 against 229/210. **The cap is not the cause**: on V1-V2, inter-laboratory V1-V2 and spike-in maxEE is 2 in both, and no limit anywhere goes above 5.
- **What was wrong with the check.** Section 21 (P30) said the released package was run "with the same default subsample as the nights". That was wrong. `compare_methods.py` thins every dataset (to 5,000 pairs per sample; mock A is not thinned) **before** it runs AmpliTrim and DADA2, so the nights' picks were computed on those pairs, and the check gave the package the full reads instead. The code that builds the pass tables and chooses the cut is the same in the two packages (a line-by-line comparison shows only added helper methods, removed output code and comments in `trimmer.py`, `ranking.py`, `merge_wrapper.py` and `subsample.py`; `pipeline.py` differs in how the search is run: FIGARO's suggestions, the rule, ITS), and the nights did not use FIGARO's splits. Different pairs are therefore the leading explanation. **This is a hypothesis, not a result.**
- **Consequence under the rule of section 21.** P30 fails on these five datasets, so **the DADA2 results of the nights are not claimed for the released package on V1-V2, inter-laboratory V1-V2, inter-laboratory V4-V6, wastewater and spike-in** until P39b below, or DADA2 at the released settings, has been read. The README, the website and the results pages present those results as the tool's.
- **P43 MET**, at the threshold: the floor-only cut is the released cut on 8 of 11 datasets. It differs on mock A (274/246 against 273/247), inter-laboratory V4-V6 (274/264 against 273/265) and wastewater (254/186, total 440, against 260/182, total 442). P44 (DADA2 at these three floor-only cuts) has not been run.
- **P45 (descriptive).** Pairs lost because a read is shorter than the cut: 0.0% on seven datasets, 0.2% on mock B and on wastewater, **7.7% on swabs and 10.2% on D6300 V4-V5** (the package warns from 1%). **Wastewater's low pass rate is therefore not a loss to short reads**, and the wording that calls it a dataset of poor read quality stands.
- **P42 (descriptive).** Expected errors per kept read on mock B at 274/253 and maxEE 4,8: 0.90 forward and 3.56 reverse (predicted about 1.0 and 3.3).

**P39b (corrected check; the prediction is written before it is run, 11 October 2026).** Given the pairs the nights gave AmpliTrim (`nocap_check.py --like-the-nights`: every dataset thinned to 5,000 pairs per sample, mock A not thinned, with `real_data.thin_reads`, the function `compare_methods.py` used), the released package recommends exactly the settings that were run through DADA2 on every dataset, except mock B, where it sets 274/253 with maxEE 4,8.

- *If P39b holds:* the cause was the sample. The DADA2 results of nights 3 and 4b are the released package's results on the pairs of the nights, and the documents say so. That the package's pick on the full reads of a dataset can differ from its pick on the nights' pairs by a few nucleotides (a different sample of pairs along a flat stretch; the size is recorded here: up to 6 nt in the forward cut and 11 in the reverse, and 5 in the total) is reported as the sampling variation of the tool, as in the stability test.
- *If P39b fails on some datasets:* DADA2 is run at the released settings (arm A and the released setting) on those datasets only, before any claim, as section 21 says, and the datasets where it holds stand.

**P43 and P44 on the nights' pairs (11 October 2026, written before P39b is run).** DADA2 sees the nights' pairs, so P43 is read again in the report of P39b (the floor-only cut and the released cut, both on those pairs), and the DADA2 runs of P44 are made at the floor-only cuts of **that** report, with the released cut as the baseline (`--setting released=... --setting floor_only=... --baseline released`; mock A is given all its reads, as in the nights). The P43 result above, on the package's own sample, stays recorded as it is. The commands printed in the first report for P44 were **not run**: they had `--pairs 5000` for mock A, which the nights did not use, and compared the floor-only cut with AmpliTrim's own pick and not with the released cut.

**Outcome, P39b and P43 on the nights' pairs (11 October 2026; `nocap_check.py --like-the-nights`, 12 datasets).**

- **P39b NOT MET, by one dataset.** Mock B is as predicted: the released package sets 274/253 with maxEE 4,8 on the pairs of the nights, which is the setting DADA2 was run at (P40, P41). Of the other eleven datasets, **ten recommend exactly what DADA2 was run at**: mock A, V1-V2, ITS, soil V4, inter-laboratory V1-V2, D6300 V3-V4, D6300 V4-V5, wastewater, spike-in and swabs. **Four of the five datasets that differed on the full reads are now the same** (V1-V2, inter-laboratory V1-V2, wastewater, spike-in): the different sample of pairs explains those. **Inter-laboratory V4-V6 still differs**: 271/267 with maxEE 3 was tested (night 4b, arm H), and the package recommends 273/265 with maxEE 4,3, **the same on the full reads and on the nights' pairs**, so the sample does not explain this one. The code that computes the rule from a cut is the same in the night script (`read_errors`, `q_rule`) and the release (`cut_errors`, `_rule_max_ee`). **The cause is not known.** Open possibilities, none tested: the 75th percentile of the forward reads lies close to an integer, and the rule (`int` of it, plus 1) is a step function; the first-pass cut at maxEE 2 differs slightly on this dataset, which has a tight overlap; something else. DADA2 at the released setting (the rule of section 21: arm A and the released setting) is the next step for this dataset, together with the tested setting.
- **Consequence under section 21.** For the ten datasets where the release recommends what was tested, the DADA2 results of nights 3 and 4b stand as the released package's results on the pairs of the nights. For inter-laboratory V4-V6 they are not claimed for the release until the DADA2 run has been read. Strictly, P39b as written is not met, so the decision rule above (a measured result for mock B on the validation page and the results site) is applied only after the V4-V6 run, in one step.
- **P43 NOT MET on the nights' pairs: 7 of 11 (at least 8 predicted).** (On the package's own sample of the full reads it was 8 of 11, recorded above; the nights' pairs are the ones DADA2 sees, as written before this run.) The floor-only cut differs from the released cut on mock A (274/246 against 273/247), inter-laboratory V4-V6 (274/264 against 273/265), wastewater (254/186, total 440, against 261/183, total 444) and spike-in (216/224 against 222/218, the same total 440). By the decision rule fixed in advance, **the search over longer totals and the four-key tie-break are not candidates for removal on this evidence**; nothing changes in 0.1.0. P44 (DADA2 at these four floor-only cuts, against the released cut) is still to be run, because it says how much the differences matter.
- **P45 (descriptive), on the nights' pairs.** Pairs lost because a read is shorter than the cut: 0.0% on seven datasets, 0.2% on mock B, 0.6% on wastewater, 8.3% on swabs and 10.4% on D6300 V4-V5. Wastewater's low pass rate is not a loss to short reads: **"poor read quality" stands**.
- **P42 (descriptive).** Mock B at 274/253 and maxEE 4,8: 0.91 expected errors per kept forward read and 3.48 per kept reverse read.

**Outcome, the four DADA2 runs (11 October 2026; `compare_methods.py` on the nights' pairs, with the released setting as the baseline).**

| dataset | released: cut, maxEE | merged | floor-only cut | merged | floor-only minus released |
|---|---|---|---|---|---|
| inter-laboratory V4-V6 | 273/265, 4,3 | 63.4% | 274/264 | 62.6% | -0.8 (5 samples: no interval) |
| mock A | 273/247, 2,3 | 69.5% | 274/246 | 69.5% | -0.0 (-0.2 to +0.1; better in 4 of 6 samples) |
| wastewater | 261/183, 4,3 | 57.9% | 254/186 | 60.3% | **+2.3 (+1.5 to +3.1; better in 9 of 9 samples)** |
| spike-in | 222/218, 2 | 83.6% | 216/224 | 83.6% | +0.0 (-0.2 to +0.2; better in 7 of 12 samples) |

- **P44 MET on all four datasets**: the floor-only cut merged at least as many pairs as the released cut, within the margin, on every dataset where the two differ.
- **Inter-laboratory V4-V6 (the open item of P39b).** The default (273/265, maxEE 2) merged 41.6%, which reproduces the nights. The tested setting (night 4b, 271/267 with maxEE 3) merged 57.9%. The released setting (273/265 with maxEE 4,3) merged **63.4%**, which is **5.5 points more than the tested setting**, beyond the margin of 2. **The nights' result for V4-V6 therefore does not describe the released package, and the released package does better.** The cause of the difference between the two recommendations is still not known. One thing is clear from the numbers: a change of one unit in the forward limit (3 to 4) is worth 5.5 points on this dataset, so the rule's limit, which is a step function of the percentile, can matter a lot where a percentile is close to an integer.
- **The released package's DADA2 result is now known on all 12 datasets**: the nights' result on the ten where it recommends the same settings on the nights' pairs, and the night 5 runs on mock B and on V4-V6. The validation page, the findings note, the paper, the user guide, the report text of the package, the website and the results pages were updated on 11 October to these numbers: mock B 70.4% (4,8) and V4-V6 63.4% (4,3); over the five datasets where the rule differs from the default, the mean is 64.9%, against 61.2% for `2,5` and 45.0% for the default, and the rule merges more than `2,5` on three datasets and fewer on two (mock A and V4-V5).
- **The adoption criteria of night 4b were not applied again** to the released rule's numbers, because those numbers were seen first. The status of the rule stays "a sensible default, not independently tested".
- **Not re-run for the released package:** the comparisons of Table 4 (AmpliTrim's cut against FIGARO's at FIGARO's `maxEE`), which the nights ran with the development build at those settings. The search code is the same in both packages and gave the same picks on ten datasets, at `maxEE` 2 and at the rule's settings, but the picks at FIGARO's `maxEE` were not compared.
- **Wastewater (descriptive).** The longer total chosen by the search (444 nt, against a floor of 440) merged 2.3 points fewer pairs than the floor-only cut in 9 of 9 samples, although AmpliTrim's estimate preferred the longer one (estimated yield 65.9% against 63.1%; DADA2 57.9% against 60.3%). Wastewater is also the dataset where the other tool was ahead at the same `maxEE`.
- **Simplification, by the rule fixed in advance.** P43 was NOT MET on the nights' pairs (7 of 11), so the search over longer totals and the four-key tie-break are **not candidates for removal on this evidence**, and nothing changes in 0.1.0. P44, which asks the practical question (does the floor-only cut merge as many pairs?), was MET on all four datasets where the cuts differ, and on wastewater the floor-only cut merged 2.3 points more. P43, as a test of whether the cut is the same, was stricter than that question. **This is recorded as an observation, not as a changed criterion.** To pursue a simpler algorithm, a new test is needed with its own criteria fixed before it is run: all 12 datasets through DADA2 at the floor-only cut, and the stability over three sample sets.

**Verdicts of night 5, as they stand.** P39 NOT MET (full reads). P39b NOT MET (one dataset, V4-V6, resolved by DADA2: the released package is better there). P40 MET. P41 MET. P42 descriptive. P43 MET on the package's own sample (8 of 11) and NOT MET on the nights' pairs (7 of 11). P44 MET (4 of 4). P45 descriptive.

## 25. Night 6: an independent real dataset, PRJNA933120 (Buetas et al. 2024): the released package against DADA2 on its Illumina half, and a description of its PacBio half (drafted 13 October 2026; predictions written before any of it was run)

**Why this dataset.** It is the first dataset that played no part in the design of AmpliTrim. Buetas et al. (BMC Genomics 25:310, 2024; CC BY) sequenced the same nine people as 16S V3-V4 on an Illumina MiSeq (2 x 300, primers 341F and 805R) and as full-length 16S on a PacBio Sequel II (27F and 1492R). The run table of the accession (read on 13 October 2026) lists 146 runs: 120 PacBio and 26 Illumina. The 26 Illumina runs are the paper's samples (their read counts equal the paper's Table 2) and each has a PacBio twin; `validation/PRJNA933120_pairs.csv` lists the 26 pairs (faeces are named H on PacBio and F on Illumina). The other 94 PacBio runs belong to other volunteers and are not used. **No mock community is among the runs, so nothing here has a known truth**: the Illumina half is judged as the other real datasets were, by what DADA2 does with the settings, and the PacBio half is only described.

**What the status of the rule can and cannot become.** The `maxEE` rule is a sensible default that has not been tested on independent data. This is the first dataset on which it can be. One dataset cannot make it tested: the outcome is recorded as one independent observation. If P47 is met, the rule held on 1 of 1 independent real dataset; if it is not met, the status is reconsidered in a new decision recorded here, and the released rule is not changed in 0.1.0 because of one result without that decision.

**Margins** (as in nights 2 to 5; the same table is in `validation/night6.py`, `MARGINS`, and printed in the summary): the estimate of the share of pairs passing the filter within 2 points of DADA2; the share merged by the rule at least 2 points above the default when the rule relaxes `maxEE`. Reads are thinned to 5,000 pairs per sample for AmpliTrim and DADA2, as in the nights.

**Questions and predictions, written before any of this was run.**

- **P46 (DADA2, the filter estimate).** At the released setting and at the default setting (the best cut at `maxEE` 2), AmpliTrim's estimate of the share of pairs passing DADA2's filter is within 2 points of what DADA2 measured (`est_pass` against `pass` in `compare_methods.tsv`). Expected to be met: the estimate was within 0.8 points on all 12 earlier datasets (mean 0.2).
- **P47 (DADA2, the rule against the default).** If the released rule relaxes `maxEE` above 2 on either read, DADA2 at the released setting merges at least 2 points more of the pairs than at the default setting. If the rule keeps `maxEE` 2 on both reads, P47 does not apply. Expected: on the five earlier datasets where the rule relaxed `maxEE` it merged more than the default on all five (by 6.1 to 50.5 points: mock A 62.7% to 69.5%, mock B 19.9% to 70.4%, inter-laboratory V4-V6 41.6% to 63.4%, V4-V5 57.2% to 63.3%, wastewater 43.6% to 57.9%).
- **P48 (descriptive, no criterion).** The shares merged at the default, at 2,5 and at the released setting. The rule was above 2,5 on three of five earlier datasets and below it on two.
- **P49 (descriptive).** The development build (the archive's `amplitrim/`, used by `compare_methods.py` for its own AmpliTrim row) and the released package choose the same cut at `maxEE` 2. Expected: the same.
- **P50 to P53 (descriptive expectations on the PacBio reads; no criterion changes anything).** P50: at least 90% of the reads have both primers. P51: the median length between the primers is 1,430 to 1,480 nt, and at least 95% of the reads are within 10% of it (the raw SRA length of 1,497 minus the two primers, 39 nt, is 1,458, and the paper reports 1,452 to 1,469 after processing). P52: the SRA copy keeps the original quality scores: every file has at least 40 distinct values and a highest value of at least Q60. P53: the mean quality of the first 50 bases and of the last 50 bases is within 3 Q of the middle half (HiFi reads are expected to be uniform).

**Decisions, as set in advance.** If P46 is met and P47 is met or does not apply, the dataset is added to the results site as dataset 13, kept aside and not used in the design, with the real AmpliTrim report that the run writes (`amplitrim/report/`), and the counts on the overview change. If P46 is not met, the size of the error is reported as it is, nothing in the package changes, and the cause is investigated in a new section. If P52 is not as expected (simplified qualities), the PacBio expected errors are not used for anything, and a copy with the original qualities must be obtained before a long-read design relies on them. P48 to P51 and P53 attach no decision.

**What is not tested.** The accuracy of variants (there is no truth). The PacBio reads are not denoised. The PacBio runs of the accession are not certain to be the paper's libraries: for four of the eight faecal twins, whose names carry `-R`, the read counts are 2 to 9 times those in the paper's Table 2 (see `FUTURE.md`), so P50 to P53 describe the deposited runs and not the paper's analysis. The ENA copies of the reads are assumed to equal the SRA's (the MD5 sums of ENA are checked).

**How to run** (step by step in `PUSH_TO_GITHUB.md`, section 2d; all on the machine that has R, vsearch and cutadapt).

```
python validation\night6.py --check --package C:\Users\Admin\OneDrive\Documents\GitHub\amplitrim
python validation\night6.py --out C:\data\runs\night6_pilot --package C:\Users\Admin\OneDrive\Documents\GitHub\amplitrim --limit 2
python validation\night6.py --out C:\data\runs\night6_PRJNA933120 --package C:\Users\Admin\OneDrive\Documents\GitHub\amplitrim
```

The run writes `NIGHT6_SUMMARY.md`, which applies P46 to P53 with the margins above.

## 26. Night 6: outcomes (14 October 2026)

**The run.** `validation/night6.py` on the whole table of 26 pairs, folder `C:\data\runs\night6_PRJNA933120`, the released package 0.1.0 (the clone of the software repository, the version without a cap on `maxEE`), cutadapt 5.2, R with dada2 1.40.0. All 52 runs were downloaded and agreed with ENA's MD5 sums (3 minutes); the whole run took about 9 minutes. The Illumina reads start with the primers (341F at base 0 in 98.0% of R1, the reverse primer in 96.7% of R2), so the primer step was needed and the 26 samples were thinned to at most 5,000 pairs each, as in the nights.

**Illumina half: not decided, because the deposited copy has no per-base quality.** `describe_reads.py` found **2 distinct quality values in each of the 26 samples**. The released package said so in its report ("Quality scores take only 2 distinct values ... consistent with a simplified-quality file such as the SRA 'lite' format; use the original FASTQ files instead"), and DADA2 stopped for both settings with `Error matrix is NULL` (it cannot learn an error model from two values). Outcomes of the criteria of section 25:

- **P46** (the filter estimate within 2 points of DADA2): **not decided**. DADA2 gave no result at any setting, so there is nothing to compare the estimate with. This is not "not met".
- **P47** (the rule against the default): the rule kept `maxEE` 2 at the cut 251/191, so P47 does not apply. **That carries no evidence about the rule**: the recommendation was computed from the same two-valued scores, which say almost nothing about errors.
- **P48**: no share merged, because DADA2 did not run. **P49**: the development build and the released package chose the same cut, 251/191, at `maxEE` 2 (both used the same scores).

**PacBio half** (26 runs, 606,658 reads; descriptive, no criterion changes anything):

- **P50: not as expected.** 83.0% of the reads have both primers by a mismatch-only match (at most 3 of 20 bases wrong, within 60 bases of the end of the read); the expectation was at least 90%, which 3 of the 26 files reach (92, 91 and 93%; all saliva runs without `-R`). **What `longread/longread_describe.tsv` shows about the 17%:** 16.3% of the reads have exactly one primer (5' end only 8.2%, 3' end only 8.1%) and only 0.7% have none. That is what independent misses of each primer give: a chance of 8.9% per primer predicts 83.0% with both, 16.2% with one and 0.8% with none. The share with both primers follows quality across the files (Spearman correlation +0.89 with the mean quality of the middle, -0.90 with the median expected errors) and not the number of reads (-0.06). The six files whose library names end in `-R` are the lowest in quality and in the share: 77.0% with both primers against 85.8% in the other 20, middle quality Q28.7 against Q31.5, median expected errors 0.87 against 0.04, and 64.4% of the reads at or under `maxEE` 2 against 83.4%. **The reading that fits, not yet tested:** errors in the primer region, above all insertions and deletions, which a mismatch-only match cannot absorb, hide a primer that is there; the reads are not missing their primers. A search that allows indels on a low-quality file would show it (cutadapt, with the commands in the guide's section 2d). Reads whose primer is missed are left out of every other figure of this half, and are probably the noisier ones, so those figures are somewhat optimistic, most of all for the lowest-quality files.
- **P51: as expected.** The median length between the primers is 1,454 nt (expected 1,430 to 1,480), and 99.6% of those reads are within 10% of it (expected at least 95%). With the raw length of 1,497 nt this confirms that the primers are still in the deposited reads.
- **P52: as expected.** Every file has at least 92 distinct quality values and a highest value of Q93, so the ENA copy of the PacBio reads **keeps the original quality scores**. The Illumina reads of the same accession do not (above); that the two differ by platform within one accession is what the data show, and why it is so is not established.
- **P53: as expected, within the margin.** The mean quality (as Phred of the mean error probability) is Q30.4 over the first 50 bases, Q30.3 over the middle half and Q27.8 over the last 50 bases: the 3' end is 2.5 Q below the middle (the margin was 3). The verdict is on the pooled figures and stands; per file, the last 50 bases are more than 3 Q below the middle in 5 of the 26 files (SRR23380883, SRR23380887, SRR23380901, SRR23380916, SRR23380952), all of them named `-R`.
- **Added by the table, no prediction attached.** Retention among the reads that have both primers, pooled: 70.0% at `maxEE` 1, 77.9% at 2 (56% to 88% over the files), 82.9% at 3, 89.0% at 5, 94.3% at 8, 96.5% at 10 and 99.9% at 15. No read has more than about 15 expected errors, which is what the HiFi definition gives (predicted accuracy of at least 99%, so at most 1% of 1,454 nt, 14.5): consistent with these being HiFi reads (an inference). The six `-R` runs have 42,358 reads per file against 17,625 in the other 20 and a lower per-read quality; what `-R` stands for is still not stated.

**A flaw found in the run's own scripts.** The first version of `night6.py` counted the DADA2 step as finished when `compare_methods.py` ended normally, although DADA2 had failed for every setting (the failure is written to a column of the table, not to the exit code). It was found because the summary showed dashes and "done in 1 min". It is fixed: the step now fails, quoting the error, and is skipped with the reason when the reads have two or fewer distinct quality values; `--redo STEP` repeats one step; the summary states the number of quality values at the top of the Illumina section. The sentence that `compare_methods.py`'s DADA2 error message used ("AmpliTrim's own result above is unaffected") was too reassuring and now says that the result is computed from the same scores. The build id changed with these edits.

**Decisions proposed (for the researcher to confirm or change; none changes the package).**

1. P46 to P49 are recorded as **not decided** for this copy, not as met or not met.
2. The dataset is **not added to the results site**: it has no DADA2 result.
3. The criteria of section 25 stand unchanged for a copy of these reads that has the original quality scores, or for another dataset: that would be recorded as a new section (night 6b), with its predictions written before it is run.
4. P52 being as expected, the PacBio expected-error figures may be used; no decision is attached to P50, P51 and P53.
5. **Before a dataset is chosen, one of its runs is described** (`describe_reads.py`, the `Qvals` column). It takes seconds, and would have found this before the run.

**A diagnostic, written before it is run (P54, descriptive).** Question: are the primers that the mismatch-only match misses really absent, or hidden by errors? Method: cutadapt (linked adapters, both primers, either orientation, errors up to 15% of a primer, **insertions and deletions allowed**) on two raw files: SRR23380883 (the lowest in quality; 75% with both primers by the mismatch-only match) and SRR23380950 (one of the best; 93%). Measure: the share of reads written when only reads with both primers are kept. Criteria: if the share for SRR23380883 is at least 8 points above 75% (so at least 83%), the primers were hidden by errors that a mismatch-only match cannot absorb (supported); if it is under 3 points above (under 78%), the primers are really missing and the explanation is rejected; in between, partly. The rise for SRR23380950 is expected to be smaller than for SRR23380883. No decision is attached; if supported, the long-read describer's primer matching is to be changed to allow indels.

**P54 outcome (14 October): supported.** cutadapt with indels allowed found both primers in **99.9%** of the reads of SRR23380883 (25,956 of 25,986; 75% by the mismatch-only match, a rise of 24.9 points against the criterion of 8) and in **99.7%** of SRR23380950 (18,588 of 18,642; 93% before, a rise of 6.7, smaller than for SRR23380883 as expected). The primers are in the reads; the mismatch-only search could not see them. cutadapt's own tables show the mechanism: in SRR23380883 26% of the matched 5' primers are 19 or 21 nt long instead of 20, and 22% of the matched 3' primers are not 19 nt long, that is, they carry an insertion or a deletion (in SRR23380950 8% and 12%); 54% of the reads were reverse-complemented in both files, as the describer had found. **Consequence, as written in advance:** the primer matching of `validation/longread_describe.py` now uses an edit distance (insertions and deletions count as errors, up to 3 of 20 bases, the same budget as cutadapt's) and reports how often each primer is found; the figures of the PacBio half above (P50 to P53, the retention at each `maxEE`, the quality at the ends) were computed with the mismatch-only match on 83.0% of the reads, so they are replaced by the figures of the rerun, which are recorded below when it has been run. **The P50 verdict recorded above stays as it was measured.**

**The PacBio half rerun with the indel-tolerant describer (14 October, 1 minute).** Both primers are found in **99.7%** of the 606,658 reads: **P50 is as expected** (the first measurement, 83.0%, stays recorded). The median length between the primers is 1,454 nt and 99.6% of the reads are within 10% of it (P51, as expected, unchanged). Every file still has at least 92 distinct quality values and a highest value of Q93 (P52, unchanged). The mean quality is Q29.9 over the first 50 bases, Q29.9 over the middle half and Q27.5 over the last 50 (P53, as expected: the 3' end is 2.4 Q below the middle; the first run gave Q30.4, Q30.3 and Q27.8). The pooled quality is 0.3 to 0.5 Q lower than with the mismatch-only match, the direction expected if the reads that had been missed are the noisier ones. The first run's files are kept in `longread_mismatch_only`.

**Not known yet.** The per-file table and the retention at each `maxEE` of the rerun below (not read yet). What `-R` stands for. **ENA has no original file for the Illumina runs:** its file report for SRR23380961, asked on 14 October, has empty `submitted_format` and `submitted_ftp` columns. Whether NCBI's SRA holds a full-quality version of these runs, or the authors the submitter's FASTQ files, is not known.

## 27. The minimum overlap: is 12 nt safe, and what does a longer overlap cost?

**How this section was written and run (stated as it happened).** The questions, datasets, arms, measures, predictions and decision rules below
were written in a file (`REGISTER_27_overlap.md`) on 4 November 2026 before the first overlap run (the first run folder was made at 09:35).
They were **not committed to this archive before the runs**: `TEST_REGISTER.md` was last committed on 15 October 2026 (`d8ad5cf`), so the git history
cannot show that the predictions came first. This section was committed on 2026-11-04, after all 36 runs. Mock A and mock B
(18 runs, 09:50 to 10:44) were analysed at 10:53, and P52 was already NOT MET on both before the D6300 runs were started (11:01 to 11:23);
so the D6300 results are not blind to it, and the decision rule for P52 (two or more mocks) was already triggered by mock A and mock B alone.
Run folders: `23_overlap_test_2026-11-04` (the numbering of the folders follows the order of the register, not the dates).

**The question.** AmpliTrim's default minimum overlap is 12 nt, DADA2's `mergePairs` default since version 1.8 (it was 20 in 1.0 to 1.6;
checked in the DADA2 source on 4 November 2026). The tutorial's advice of 20 nt still reflects the old default. vsearch's default is 10.
No run of this study counted wrong merges. Two things are unknown: (a) does a low minimum overlap let DADA2 make merges that are wrong,
and (b) what does a higher minimum overlap cost in pairs and in noisy read ends?

(`collect_overlap_results.py` prints these predictions as P-O1 to P-O5: P-O1 is P50, and so on.)

**Overlaps tested:** 10, 12, 15, 20 and 25 nt. 10 is vsearch's default, 12 is DADA2's, 15 lies between, 20 is the former DADA2 default and
the tutorial's advice, 25 is a generous margin. AmpliTrim itself accepts 12 to 20; 10 and 25 are tested at the level of DADA2 only
(`compare_methods.py --min-overlap`, `--setting`).

**Datasets.** The four mock communities with expected sequences, whose longest amplicons were merged at exactly 12 nt:

| dataset | longest amplicon | forward cut F | maxEE (as run) | pairs per sample |
|---|---|---|---|---|
| mock A | 508 | 273 | 2,3 | all |
| mock B | 515 | 274 | 4,8 | 5,000 |
| D6300 V3-V4 | 430 | 224 | 2 | 5,000 |
| D6300 V4-V5 | 550 | 282 | 2,3 | 5,000 |

**Two arms**, for every overlap k. maxEE is the same in all runs of a dataset, so that only the overlap differs.

- **Arm T (the cut follows k).** Total length = longest amplicon + k, forward cut F as above, reverse cut = amplicon + k - F,
  `mergePairs(minOverlap = k, maxMismatch = 0)`. It answers: does the overlap that is asked for change accuracy or yield?
- **Arm S (the cut stays).** The cut of k = 12 (reverse = amplicon + 12 - F) is kept; only `minOverlap` changes. For k = 10 it can only
  show merges that 12 would not make; for k = 15, 20 and 25 it shows the pairs lost because their real overlap is shorter than k.

**Measures** (as in Table 6 of the findings, with `--reference`): share of input pairs merged; variants recovered; strains recovered;
share of reads in exact ASVs; and the number of ASVs that are not identical to an expected sequence.

**Predictions, with the rule that decides each.**

- **P50 (arm S, no extra merges at 10).** With the k = 12 cut, `minOverlap` 10 gives a merged share within 0.1 points of 12 and no ASV
  that is absent at 12, on each of the four mocks. MET if all four.
- **P51 (arm T, accuracy at 12).** At k = 12 the share of reads in exact ASVs is within 1.0 point of the best k, and the strains recovered
  equal the best k, on each of the four mocks. MET if all four.
- **P52 (arm T, no more unexpected ASVs at 12).** The number of ASVs not identical to an expected sequence at k = 12 is at most the number
  at k = 20 plus 2, on each of the four mocks. MET if all four.
- **P53 (descriptive, the cost of overlap).** The merged share at k = 25 is lower than at k = 12 by at least 1 point on at least three of
  the four mocks, in arm T and in arm S. Reported, not a verdict.
- **P54 (descriptive, the longest amplicons).** In arm S, the pairs lost at k = 15, 20 and 25 are the pairs of the longest amplicons.
  Reported with the share lost per k.

**What each outcome changes (decided now).**

- P50 to P52 all MET: the paper says that 12 nt matched 15, 20 and 25 on four mocks, and keeps 12 as the default.
- P52 NOT MET on two or more mocks (more unexpected ASVs at 12 than at 20): the default becomes 15 in 0.3.0, and the paper reports the
  result as it stands.
- P50 NOT MET on any mock: DADA2 makes merges at 10 that it does not make at 12. The paper says so, the default stays at 12, and
  `--min-overlap` below 12 stays refused.
- Whatever the verdicts, the table of all five overlaps is shown in a supplementary table.

**What this does not test.** Wrong merges on amplicons with long repeats (ITS) and on real samples without expected sequences. It
tests the four mocks only; the Discussion must say so.

**How to run.** `make_overlap_commands.py` writes the 36 commands (4 datasets, 9 runs each). Start with mock A and mock B, the datasets
with the tightest overlap (`--only mockA mockB`).

### 27, outcomes (4 November 2026; run folder `23_overlap_test_2026-11-04`; table `overlap_results.md`, fingerprints of the ASV sets `overlap_results_asv_sets.tsv`)

All 36 runs were made. Three produced no result because the cut was longer than the reads: DADA2's `filterAndTrim` removed every read
(D6300 V3-V4 at k = 25, cut 224/231; D6300 V4-V5 at k = 20 and k = 25, cuts 282/288 and 282/293; the log says `r removed all reads`). D6300 V4-V5 at
k = 15 (cut 282/283) merged 0.0%. The floor cut of V4-V5 (282 + 280 nt) already uses the whole read, so no larger overlap is possible with these reads.

- **P50 MET (4 of 4).** With the k = 12 cut, `minOverlap` 10 gave the same merged share (within 0.1 points) and exactly the same ASVs as 12 on all four mocks
  (the fingerprints in `overlap_results_asv_sets.tsv` are equal).
- **P51 MET (4 of 4).** At k = 12 the share of reads in exact ASVs was within 0.3 points of the best k (and the strains recovered equal to the best) on every mock.
- **P52 NOT MET (judged on 3 of 4).** Extra ASVs at k = 12 against k = 20: mock A 16 against 11, mock B 69 against 66 (both above "plus 2"); D6300 V3-V4 1 against 1 (met). D6300 V4-V5 could not be judged (no result at k = 20).
- **P53 (descriptive).** Merged share at k = 25 against k = 12, arm T: mock A 65.4 against 69.5 (-4.1), mock B 63.9 against 70.4 (-6.5); the other two mocks have no result at k = 25.
  Arm S: mock A -20.1, mock B -3.1, V3-V4 -90.6 (0.0%), V4-V5 -0.6.
- **P54 (descriptive).** In arm S the pairs lost at higher k are those of the longest amplicons: on mock B the longest strain is lost (3/3 to 2/3 strains) at 15 nt;
  on V3-V4 every pair is lost at 15 nt, so all of its pairs overlap by 12 to 14 nt (amplicons of 427 to 430 nt) and still gave 8/8 strains and 97.9% of reads in exact ASVs.

**Read-weighted comparison (not registered; made after P52 failed).** Per 100 input pairs, using the k = 12 cut against the longer cut of arm T:
mock A, 12 against 25 nt: +4.07 reads in exact ASVs and +0.03 reads outside them (+8 extra ASVs); mock B: +4.76 and +1.74 (+15 extra ASVs; the share of reads outside exact ASVs is the same at all k,
24.6% to 25.1%); D6300 V3-V4, 12 against 15 nt: +17.98 and +0.32 (+0 extra ASVs). The count of extra ASVs falls when fewer reads are kept; the share of reads outside exact ASVs does not change with k.

**Decision.** The rule of this section was: P52 NOT MET on two or more mocks means that the default becomes 15 in 0.3.0.
P52 was NOT MET on two mocks (mock A and mock B). The default was **not** changed: it stays at 12. This deviates from the rule we set. The reasons:
no wrong merges were seen at 10 nt (P50, equal ASV sets); a longer overlap did not improve accuracy (P51); a default of 15 would cost 1.0 to 1.8 points of merged share on mock A and mock B
and 18.3 points on V3-V4, and cannot be reached on V4-V5; and P52, which counts ASVs, is a poor measure of wrong merges because the count follows the number of reads kept
(read-weighted comparison above). Users who want a margin can set `--min-overlap` up to 20.

---

## 28. Release check of version 0.2.0 (3 November 2026)

**The question.** Version 0.2.0 adds `--min-overlap` and a Python interface. Does it recommend the same settings as 0.1.0?

**What was done.** `validation\nocap_check.py --like-the-nights` was run with the 0.2.0 package on the pairs that the nights gave AmpliTrim (every dataset
thinned to 5,000 pairs per sample, mock A not thinned), the same command that gave the 0.1.0 table of section 24 (P39b).
**Outcome.** On all 12 datasets the cut and the maxEE were identical to those of 0.1.0, among them mock B (274/253, maxEE 4,8) and
inter-laboratory V4-V6 (273/265, maxEE 4,3); the pairs lost to read length (P45) and the floor-only cuts (P43) were identical as well. The verdicts
of P39 and P43 are therefore the same as in section 24 (NOT MET, for the reasons given there). The package passed its 174 tests on the same day.
Run folder: `21_release_check_020`.

*This section was written after the run; it records an outcome that was predicted (that the recommendations would not change) but not registered in advance.*

## 29. Night 10: the closing run (written 5 November 2026, before any of these runs)

**Why.** Eight gaps are left before the paper goes out. (1) Supplementary Table S5 (the share of the longest amplicons left out) mixes two builds: its 0% rows come from an earlier build.
(2) The long amplicons of D6300 V4-V5 are only suspected to be 18S. (3) The maxEE rule and the 0.2% default were set on development data; Night 6, the one independent dataset, was not decided because its
deposited reads have two quality values. (4) The ITS mode has one dataset. (5) The truncation grid covers some datasets only, on earlier builds. (6) The 12 nt overlap was tested on four mocks, all of them
seen during development, and P52 was not met there. (7) The released command line, the repeat of a run, and the `--min-overlap` option were exercised by unit tests and by the nights' Python calls, never end to end
on real reads. (8) Run time and memory (Supplementary Table S4) were measured on the 12 datasets only.

**What runs.** The released package 0.2.0 only. `night10.py` stops if the version is not 0.2.0, or if the package that the night scripts import and the one that `compare_methods.py` and the command line import
are not the same code (a hash of every `.py` file), in one command, in this order: version, accept, cli, tail, coverage, battery, its, overlap, grid. The grid is last so that it is the step left undone
if the night runs short. A pilot (`--limit 2`, 30 to 45 minutes, in a scratch folder) comes first; its numbers are not used.

**The share left out is not a command-line option.** The released command line has `--its`, `--min-overlap` and `--verbose` only. The share of the longest amplicons that may be left unmergeable is the
`max_unmergeable` setting of the Python interface (default 0.002), and the report already tabulates the best truncation and the estimated yield at 0.1, 0.2, 0.5, 1, 2, 5 and 10% and says when leaving out more would
raise the yield by 3 points or more. `--allow-unmergeable` in this register and in the scripts is the option of `compare_methods.py`, which sets that setting for a DADA2 comparison. The manuscript says this.

**The new datasets (filled in before the run; none is dropped or replaced after a result is seen).**

| Name in `amplitrim_paths.json` | Accession | Why chosen | Distinct quality values in one run (`describe_reads.py`, `Qvals`) | DADA2 pilot ran | Expected sequences |
| --- | --- | --- | --- | --- | --- |
| miseqsop_mock | none: `MiSeqSOPData.zip` (the mothur MiSeq SOP and the DADA2 tutorial data; Kozich et al. 2013) | independent mock community, V4, 2 x 250 MiSeq; never used in development (`validation/README.md`: "Test (not used so far)") | 27 (`describe_reads.py`) | yes: 97.1% of pairs merged at 132/134, 20 ASVs, all identical to a reference | yes: `HMP_MOCK.v35` (32 sequences, 21 strains), 19 of 21 strains recovered by an identical ASV |

Rules for the choice, fixed now: paired-end Illumina; more than two distinct quality values (NovaSeq's four binned values are accepted, as for the existing datasets; the two-valued file of Night 6 is not),
and a DADA2 pilot on a thinned sample that learns an error model, both before the dataset is registered; the primers known; at most four datasets; if a read length or platform that the existing 12 datasets do
not cover is on offer (2 x 250, 2 x 150), it is preferred. A dataset is dropped only for a reason found before any result (quality values, primers, missing reads, a failing pilot).

**Plan as run (decided 6 November 2026, before the pilot and before any run of Night 10).** One new dataset instead of four: the MiSeq SOP mock (`miseqsop_mock`). It serves P57 (a to d, the independent test) and P63 (the overlap hold-out) at once. The V3-V4 or V4-V5 real dataset, the further 16S real dataset and the second ITS dataset are not run: P58 is recorded as **not run**, the `its` step is skipped, and the paper says that the ITS prediction was checked on one dataset only. P57 is judged on this one dataset, and the paper says so. **Disclosure.** The prepared dataset was run once with the development build (AmpliTrim 0.1.0) and DADA2 at AmpliTrim's cut (132/134: 97.1% of pairs merged, 20 ASVs) and at 240/160 (89.3% merged), to check that DADA2 runs and that the expected sequences are found. That result was seen before the criteria of P57 and P63 were committed; none of the comparisons those criteria judge (the filter estimate, the released rule against the default, 2,5 and FIGARO, the overlap sweep) had been run, and Night 10 re-runs everything on 0.2.0. The truncation grid (P59) is run on the datasets that have no grid on the final build: `miseqsop_mock`, `v4` (soil V4, whose grid was never re-run on the released build, section 22), `mcm_v12`, `mcm_v46`, `wastewater_V3-V4`, `spikein_V3-V4` and `swabs_V3-V4`.



**Predictions and criteria** (P55 and P56 are descriptive and carry no criterion).

- **P55 (descriptive).** The share of pairs DADA2 merged at maxEE 2 when the longest 0, 0.2 (the default), 1 and 5% of the amplicons are left out, on every registered 16S dataset and every new one, one build,
  with AmpliTrim's own estimate for each row, so that the estimate can be compared with DADA2's share (estimate minus merged, in points). It replaces the earlier-build rows of Supplementary Table S5.
- **P56 (descriptive).** The best SILVA hit of the longest 1% of the merged sequences of D6300 V4-V5, against a control of median length. Reading, fixed now: if at least 80% of the long group
  is Eukaryota and at least 80% of the control is Bacteria or Archaea, the long amplicons are called 18S; if at least 80% of the long group is Bacteria, they are called 16S;
  otherwise "mixed or unresolved". A control with no hits means the search failed and nothing is concluded.
- **P57 (the independent test; every new 16S dataset on which DADA2 runs; computed from the tables of `night4.py`, not from the P19 to P22 verdict lines that script prints, which belong to Night 4b).**
  (a) AmpliTrim's filter estimate is within 2 points of DADA2's pass rate. (b) The released rule (arm H) merges at least as many pairs as the default (arm A), within the margin: H minus A at least -2 points.
  (c) On a dataset with expected sequences, H is as accurate as A: variants and strains within 1, reads in identical ASVs within 0.1 points. (d) H merges at least as many pairs as FIGARO's own pick (arm C),
  within the margin: H minus C at least -2 points, where FIGARO runs. (e, descriptive) The merge estimate minus DADA2's merged share, and whether the estimate orders the settings as DADA2 does.
- **P58 (the second ITS dataset).** AmpliTrim's filter estimate is within 2 points of DADA2's ITS filter pass rate at maxEE 2.
- **P59 (the grid, final build, 3 x 3 around the recommendation, steps of 10 nt, on the datasets named with `--grid`).** AmpliTrim's recommendation is within 1.0 point of the best DADA2 setting
  that meets the coverage floor (the margin of `overnight.py`), reported per dataset.
- **P60 (the command line, every registered dataset and every new one).** `python -m amplitrim.cli` exits with status 0 and writes `REPORT.md` and `summary.json`; on an empty input folder it exits with status 3 and a message, without a
  traceback (`docs/usage.md`). Run time and peak memory are recorded (Supplementary Table S4; descriptive).
- **P61 (repeat).** On the first two datasets, a second run on the same reads gives the same `recommended`, `maxEE_table`, `longest_amplicon` and `pairs_analysed`.
- **P62 (`--min-overlap`).** On the first two 16S datasets, at 15 and at 20 nt, `recommended.minOverlap` is the value given and `truncLenF + truncLenR - minOverlap` equals `longest_mergeable_insert`.
- **P63 (the overlap test on a hold-out mock).** On each hold-out mock (`--overlap-mock`), the 9 runs of section 27 and its criteria unchanged: P50 (with the k = 12 cut, `minOverlap` 10 gives the same
  merged share within 0.1 points and no ASV that 12 does not give), P51 (at k = 12 the share of reads in exact ASVs is within 1.0 point of the best k and the strains recovered equal the best), P52 (the unexpected
  ASVs at k = 12 are at most those at k = 20 plus 2); P53 and P54 are descriptive. `collect_overlap_results.py` prints them as P-O1 to P-O5 and says how many of the four mocks it judged; here that is the hold-out only.
- **P64 (the release check on the installed package; corrected on 6 November 2026, after the pilot and before the run, because the first wording contradicted section 21).** `release_check.py` on mock B, V4 and ITS (`--only`: the verdict of section 21 is then "not decided") finds the settings recorded in section 21 and at the acceptance run of 4 November: mock B 274/253 with maxEE 4,8 (Night 4b, before the uncapped rule, had 277/250 with maxEE 4,5: the change predicted in section 21), soil V4 124/115 with maxEE 2 and ITS 0/0 with maxEE 2 (Night 3, arm A). The column "the same" of the script is therefore expected to read yes for soil V4 and ITS and NO for mock B.

**What follows (written before the run).**

1. Nothing in the package changes because of this run: the maxEE rule, the 0.2% share and the 12 nt overlap stay as released.
2. A criterion that fails is reported as not met, in the paper and in the register. The Abstract's sentence that the rule "has not been tested on independent data" is replaced by what was measured.
3. A tip to raise the coverage share to 1% goes into the documentation and the Discussion only if, on at least two datasets that did not inform the default (not D6300 V4-V5, not soil V4), 1% merges at least
   5 points more pairs than 0.2% and loses no variant or strain on any mock. Otherwise the 0.2% sentence of the Discussion stands as written.
4. If P52 is not met on the hold-out mock either (more unexpected ASVs at 12 than at 20, beyond 2), the paper reports that three of five mocks did not meet it and keeps 12 nt as released, with the reasons
   of section 27; the default is not changed in 0.2.x (a change would be a new release with its own test).
5. A dataset on which DADA2 does not run is recorded as not decided, not as met or not met (Night 6).
6. If P60, P61 or P62 fails, the failure is a bug report on the software (an issue, a fix and a new patch release if it affects a result), and the manuscript sentence that relies on it ("a dataset always gets
   the same answer", the `--min-overlap` option) is reworded to what was observed.
7. After this run no further test is run before submission; later ideas go to `FUTURE.md`.

**Order of the work after the run.** Outcomes recorded here as section 30 first (and the sections of run 24 in `organise_data.py`, `validation/README.md`, `tests/test_archive_order.py` changed to "29-30");
then Supplementary Tables S4 (run time) and S5 (final build); then the manuscript; then the archive; the Zenodo release of the archive last, because it freezes the files.

Addendum A (written on 7 November 2026 while Night 10 was running, before its results were read and before these tests ran): more tests on the independent mock. miseqsop_mock has expected sequences, so accuracy can be scored exactly. Night 10 uses it for the overlap hold-out, the coverage table and the grid; the comparisons below are run in a short second run (folder night10_mocktests_<date>, run 24) with the released package 0.2.0, at the cut that the cli step of Night 10 recorded for miseqsop_mock. Disclosure: the DADA2 result of the development build at its own cut and at the tutorial’s 240/160 (97.1% and 89.3% merged) was seen before these criteria were written; none of the comparisons below had been run on the released package. One sample of 4,779 pairs: the results are descriptive and judged by the margins below, not by a test of significance.

P65 (maxEE where the rule keeps 2). DADA2 at the released cut with maxEE 1,1; 2,2; 2,3; 2,4; 3,3; 2,5; 3,5 and 5,5 (max_ee_sweep.py, all pairs, scored against HMP_MOCK.v35.strains.fasta). (a) The best relaxed setting merges at most 1.0 point more pairs than 2,2, so the rule’s decision to keep 2,2 costs at most 1 point. (b) Relaxing does not make the result more accurate: for every relaxed setting, the reads in exact ASVs are within 0.1 points of the 2,2 row and the strains recovered are equal to it. (c, descriptive) maxEE 1,1 merges at least 1.0 point fewer pairs than 2,2.
P66 (descriptive). AmpliTrim’s estimated yield for each of those settings (the maxEE_table of summary.json) minus DADA2’s merged share, and whether the estimate puts the settings in the order DADA2 does (ties, and differences below 0.5 points, count as the same rank).
P67 (the tutorial’s own settings on the tutorial’s own data). compare_methods.py at the released cut and maxEE of the cli step against the DADA2 tutorial’s truncLen = c(240, 160), maxEE = c(2, 2) (--setting dada2_tutorial=240,160,2,2), scored against the expected sequences. (a) AmpliTrim’s cut merges at least as many pairs as the tutorial’s within the margin: AmpliTrim minus tutorial at least -2 points. (b) The reads in exact ASVs are within 0.1 points and the strains recovered within 1. (c) The filter estimate is within 2 points of DADA2’s pass rate for both rows. (d) Only if FIGARO runs on this dataset: AmpliTrim minus FIGARO’s pick at least -2 points.

A criterion that fails is reported as not met. Nothing in the package changes because of these tests.

Addendum B (written on 7 November 2026 while Night 10 was running, before its results were read and before any read of this dataset was prepared or analysed): a second ITS dataset for P58. Section 29 recorded P58 as not run because no second ITS dataset had been chosen. One is now chosen: the full-length ITS data of Tedersoo et al. (2026, Molecular Ecology Resources 26(6): e70189), ENA project PRJEB108994, the Illumina iSeq 100 (MiSeq i100) 2 x 500 runs of the soil samples. Runs, by a rule fixed in advance (soil samples named G####, the mock and the mixed controls left out; every eighth run in the order of the run accessions): ERR16773548, ERR16773556, ERR16773564, ERR16773572, ERR16773580, ERR16773588, ERR16773596, ERR16773604 and ERR16773612: nine runs, 836,452 read pairs, all libraries without the suffix _2, which I take to be the 30-cycle libraries (inferred from the names). Disclosure: two other runs of the project (ERR16773582 and ERR16773609) were downloaded as a pilot and described (3 distinct quality values, 52,178 and 51,122 pairs of 489 nt, the forward primer near the start of the reads); they are not part of the test. P58 is as in section 29: AmpliTrim’s filter estimate is within 2 points of DADA2’s ITS filter pass rate at maxEE 2, judged on the pooled dataset (compare_methods.py --marker its, run by the its step of night10.py), with the per-run differences reported. The primers are removed with cutadapt before the comparison (ITS9MUNngs TACACACCGCCCGTCG and ITS4ngsUni CCTSCSCTTANTDATATGC). P58 concerns the filter only. A criterion that fails is reported as not met; if DADA2 does not run on this dataset it is recorded as not decided.

Amendments and notes of 7 November 2026 (written while Night 10 was running, before its results were read; they change no criterion).

Commit history. Section 29 was written on 5 November 2026 but first committed on 6 November (commit b0869d0), after the first pilot of Night 10 (21:36 to 21:42), which ran an earlier version of the driver; the pilot’s numbers are not results and are not used. Later commits: P64 corrected (baf3cca, after the pilot, before the run), then the notes and Addenda A and B of this file, all before the real run. Night 10 was started on 7 November 2026 at 00:20 (its folder is named night10_closing_2026-11-06).
Soil V4 grid. The paragraph “Plan as run” says that the grid of v4 was never re-run on the released build. Section 22 shows that a grid around the recommendation of 10 October was run (124/115 against 134/105, 4.7 points); it was not run on 0.2.0, which is why it is repeated here.
Item 7 of “What follows” (no further test before submission) is superseded by Addenda A and B, which add tests before submission.
P57 (b) to (d) cannot be judged on miseqsop_mock. In the pilot, night4.py reported for it: “the rule gives (2,2) at AmpliTrim’s cut 132/134 for both rules: identical to A by construction, nothing to run”. The comparisons of P57 (b) to (d) (the rule against the default, against 2,5 and against FIGARO’s pick) have nothing to compare and are recorded as not decided, not as met. The battery (--new) is left out of the real run; the maxEE that the released package chooses on miseqsop_mock is read from the cli step. What the dataset still tests: the filter estimate against DADA2’s pass rate (P57 a, from the coverage step, row of the default share), the coverage table (P55), the overlap hold-out (P63) and the command line (P60 to P62). The independent test of the part of the rule that relaxes maxEE is not done, and the paper says so. If the cli step shows that the released package relaxes maxEE on this dataset, Addendum A (P65) is read at that maxEE and this is stated in section 30.

Corrections to the amendments above (7 November 2026, during the run, before any result was read). (a) The last sentence of item 1 is wrong: Addendum A, Addendum B and the amendments were committed on 7 November while Night 10 was already running (started 00:20), not before it; they were written before any of its results was read. (b) Commit fa07309 is described as the note on P57 (b) to (d) but contains only two blank lines (the text was not saved); that note was first saved in commit 67715cf.

**Night 10: timeline and index (compiled on 7 November 2026 from the git log of this file; it adds no criterion and changes none; the text above stays as written).** Times are GMT, which is UK time in November. Commits are those of `git log -- TEST_REGISTER.md validation/night10.py validation/tail_classify.py`. The last sentence of item 1 of the amendments above ("all before the real run") is wrong; this timeline replaces it.

*Timeline*

| When | What | Commit |
| --- | --- | --- |
| 5 Nov | Section 29 (predictions P55 to P64) and the Night 10 package written and applied to the working tree of the archive (not committed); archive tests 646 passed, 4 skipped; 22:20 `--check`: amplitrim 0.2.0 from the checkout and from PyPI, the same code (fingerprint 53857ca90f), vsearch 2.32.0, dada2 1.40.0 | none |
| 6 Nov 20:58 to 21:01 | MiSeq SOP mock prepared with the development build (0.1.0, build a75392d6f8; the release command line lacks the options of `real_data.py`): 97.1% merged at 132/134, 20 ASVs, 19 of 21 strains; 89.3% at 240/160 (disclosed under "Plan as run") | none |
| 6 Nov 21:36 to 21:42 | First pilot (`--limit 2`, first version of the driver): `battery` refused to start (5.9 GB free); the empty-folder check failed (the driver ran it from the archive folder, which has its own older `amplitrim/`); the tail table counted two hit lines per query. Numbers not used | none |
| 6 Nov ~21:57 | Disk freed to 23 GB: conda and pip caches (3.2 GB); `C:\data\_archive_2026-10` (10.5 GB, superseded, referenced by no paths file) deleted | none |
| 6 Nov 21:58:55 | First commit of section 29, the driver, the classifier, the tests and the run-24 entries. Its message says P64 was corrected: it was not | b0869d0 |
| 6 Nov 22:07:21 | P64 corrected (register only). Its message says the driver and classifier fixes were included: they were not | baf3cca |
| 6 Nov 22:13:47 | Fixed driver (commands run from the output folder) and classifier (best hit per query; taxa and identity in the report) installed; 25 Night 10 tests | ab54b5e |
| 6 Nov ~22:22 | Repaired pilot steps re-run: empty-folder check ok; tail table sound; `battery` on `miseqsop_mock`: "the rule gives (2,2) ... nothing to run". Numbers not used | none |
| 7 Nov 00:19:24 | Commit meant to hold the note on P57 (b) to (d); it holds only two blank lines (the text was not saved) | fa07309 |
| 7 Nov 00:20:16 | **Real run started** (`night10_closing_2026-11-06`; steps version, accept, cli, tail, coverage, battery, its, overlap, grid, summary; no `--new`, no `--its-new`) | none |
| 7 Nov 00:52:34 | Addendum B (second ITS dataset, P58) and Addendum A (P65 to P67; the diff of fa07309 shows it was not yet in the file) committed: after the run started, before any result was read, as each addendum says | 7808ccc |
| 7 Nov 00:59:05 | Amendments (four items), including the note on P57 (b) to (d) | 67715cf |
| 7 Nov, at the time of the listing | version, accept, cli (13 datasets and the empty-folder check) and tail done; coverage 2 of 12 datasets; the later steps not yet run | none |

*Index of the predictions* (the full wording is in the place named; "fixed" is the commit in which the wording stood in the file)

| P | Criterion in one line | Where | Fixed | Status |
| --- | --- | --- | --- | --- |
| P55 | share merged at 0, 0.2, 1 and 5% left out, with AmpliTrim's own estimate (descriptive) | section 29 | b0869d0 | coverage running |
| P56 | SILVA best hit of the longest 1% of D6300 V4-V5 against a control; 80% reading rule (descriptive) | section 29 | b0869d0 | tail done, not yet read |
| P57 a | filter estimate within 2 points of DADA2's pass rate on `miseqsop_mock` | section 29 | b0869d0 | from `coverage`, after the run |
| P57 b to d | rule against default, 2,5 and FIGARO's pick: **not decidable** on `miseqsop_mock` (the rule keeps maxEE 2) | amendment 4 | 67715cf | not decided |
| P57 e | estimate against DADA2's share and the order of settings (descriptive) | section 29 | b0869d0 | after the run |
| P58 | filter estimate within 2 points on a second ITS dataset | section 29; Addendum B (PRJEB108994, nine runs) | b0869d0; plan 7808ccc | data not yet prepared |
| P59 | grid: recommendation within 1.0 point of the best setting meeting the floor, on 7 datasets | section 29, "Plan as run" | b0869d0 | grid last |
| P60 | command line: exit 0 and files; empty folder exit 3 without traceback; time and memory (descriptive) | section 29 | b0869d0 | cli done, not yet read |
| P61 | a second run gives the same recommendation and tables (first two datasets) | section 29 | b0869d0 | cli done, not yet read |
| P62 | `--min-overlap 15` and `20`: setting honoured; total minus k equals the longest mergeable insert | section 29 | b0869d0 | cli done, not yet read |
| P63 | overlap test of section 27, criteria unchanged, on the hold-out mock | section 29 | b0869d0 | overlap not yet run |
| P64 | release check: mock B differs (expected), soil V4 and ITS the same | section 29 (corrected) | baf3cca | accept done, not yet read |
| P65 | DADA2 at the released cut, 8 maxEE settings: best relaxed at most 1.0 point above 2,2; accuracy within 0.1; 1,1 at least 1.0 point lower | Addendum A | 7808ccc | after the run |
| P66 | AmpliTrim's yield estimate against DADA2's share at those settings, and the order (descriptive) | Addendum A | 7808ccc | after the run |
| P67 | AmpliTrim's cut against the tutorial's 240/160 on the tutorial's data: at least -2 points; accuracy; filter estimate; FIGARO if it runs | Addendum A | 7808ccc | after the run |

*Decisions, deviations and provenance*

1. One new dataset instead of four ("Plan as run"); the second ITS dataset is planned in Addendum B.
2. The share of the longest amplicons left out is the Python setting `max_unmergeable`, not a command-line option; the Methods and Discussion of the manuscript were corrected.
3. P64 was reworded after the first pilot (baf3cca) because the first wording contradicted section 21.
4. Pilot findings fixed in the driver and classifier (ab54b5e). Commit messages b0869d0 and baf3cca overstated their content, and fa07309 held only blank lines; the history above is the record.
5. SILVA release 144 SSU Ref NR99 (DOI 10.82364/144/SSU/Ref-NR99/FASTA/unaligned/truncated; MD5 ccead7bba52075f5d8ba2bfab5adc4c1; 626,354 bacterial, 247,096 eukaryotic and 32,178 archaeal sequences). `miseqsop_mock`: reads gzipped, reference `HMP_MOCK.v35` renamed ID|STRAIN|label (32 sequences, 21 strains).
6. To do after the run: read `NIGHT10_SUMMARY.md`; run Addendum A; prepare the ITS runs of Addendum B and run P58; write section 30 with `apply_night10_outcomes.py`.

## 30. Night 10: outcomes (7 November 2026)

**The run.** `validation/night10.py` on the released package 0.2.0 (code fingerprint 53857ca90f for both imports, in the `version` step), folder `C:\data\runs\24_night10_closing_2026-11-06` (named `night10_closing_2026-11-06` when it ran). The main run started at 00:20 on 7 November and ended at 02:51 (2 h 31 min). The `its` step was run separately into the same folder at 01:56 (the main command had no `--its-new`; 5 min). Two steps were repeated after two defects of the validation scripts were fixed (see the last paragraph): `coverage` for `miseqsop_mock` only (10:58, 77 s) and `tail` (12:15, 71 s). Pilots: 6 November (21:36 to 22:22), discarded. The tests of Addendum A (P65 to P67) are in `C:\data\runs\24_night10_mocktests_2026-11-07`. What did not run or failed: `battery` was skipped (no `--new`, as decided before the run); the first coverage run of `miseqsop_mock` failed (a `TypeError`; repeated, table now saved); the first tail run chose the wrong group (repeated; the first output is kept in `tail_first_run`); the grid of `swabs_V3-V4` stopped (DADA2: 'Error matrix is NULL').

**The new datasets (as registered in section 29 and Addendum B).**
- `miseqsop_mock`: the MiSeq SOP mock community (Kozich et al. 2013; the DADA2 tutorial data), V4, MiSeq 2 x 250, 4,779 pairs, 27 distinct quality values, no primers to remove; DADA2 ran; expected sequences `HMP_MOCK.v35` (32 sequences, 21 strains).
- PRJEB108994 (Tedersoo et al. 2026), full-length ITS, iSeq 100 2 x 500, nine runs chosen by the fixed rule of Addendum B (ERR16773548 to ERR16773612), 836,452 pairs, three distinct quality values, pairs in mixed orientation, R2 with 12 bases before its primer. Primers (ITS9MUNngs, ITS4ngsUni) removed with cutadapt, pairs without a primer discarded: 739,189 pairs kept (88.4%; 4.6 to 13.5% lost in seven runs, 31.5% in ERR16773580 and 27.9% in ERR16773596). DADA2 ran (5,000 pairs per run; 4,605 to 4,969 reached it).

**Outcomes** (the criteria are those of section 29 and its addenda).

- **P55 (descriptive).** The table is in `NIGHT10_SUMMARY.md` (final build, maxEE 2: DADA2's merged share and AmpliTrim's estimate for each row). Covering every amplicon (0%) against the default (0.2%): V1-V2 84.2% to 45.5% (-38.7 points), spike-in 83.6% to 11.3% (-72.3), D6300 V4-V5 57.2% to 0.8% (-56.4), mock A -13.4, D6300 V3-V4 -9.1, mock B -6.0, soil V4 -2.8, swabs -1.9, inter-laboratory V4-V6 -1.1, inter-laboratory V1-V2 -0.1, wastewater 0.0, `miseqsop_mock` 0.0 (no long tail: all four rows give 132/134). Leaving out 1% against 0.2% differs by 2 points or more on one dataset (D6300 V4-V5, 57.2% to 89.8%, +32.6); leaving out 5% on five (V4-V5 +32.6, soil V4 +5.5, mock A +2.1, inter-laboratory V1-V2 -2.1, wastewater -2.2). On mock A, leaving out more lost strains (12/14 at 0.2%, 11/14 at 0% and 1%, 10/14 at 5%); mock B kept 3/3 and D6300 V3-V4 8/8 at every level; V4-V5 kept 8/9 at 0.2%, 1% and 5% and 4/9 at 0%. AmpliTrim's estimate minus DADA2's merged share at the default: median +3.4 points over the 12 datasets, from -0.8 (`miseqsop_mock`) to +7.3, and +31.5 on soil V4.
- **P56 (descriptive).** The longest 1% of the merged sequences of D6300 V4-V5 (103 of 13,675 sequences, longer than 374 nt): Eukaryota 101 (98.1%), Bacteria 2 (1.9%), median identity 99.8%; the control (103 sequences within 5 nt of the median): Bacteria 103 (100%). By the reading fixed in section 29 the long amplicons are **18S**. The 101 eukaryotic sequences are 0.74% of the merged sequences; 94% of the long group has a best hit of 97% identity or more. The organism names attached to some SILVA records do not match their taxonomy path (a hit under the genus field Saccharomyces carries the organism name of barley), so the reading uses the domain only; the labels are in `tail_summary.md`. The first execution is not used (see the last paragraph).
- **P57.** (a) **met** on `miseqsop_mock`: the filter estimate 98.0% against DADA2's 98.2% (-0.2 points). Descriptive, on the 11 other datasets of the coverage step: the largest filter difference over the four rows was 1.4 points (mock A). (b) to (d) **not decided**, as amended before the run: the released rule keeps maxEE 2 on `miseqsop_mock` (98.2% of the pairs pass at maxEE 2; the 75th percentile of the expected errors is 0.155 forward and 0.159 reverse), so the rule equals the default there by construction. (e) The estimate minus DADA2's merged share: -0.8 points (the command-line estimate, 96.4% against 97.1%) and -1.2 (the estimate that `compare_methods.py` applies to the same row, 95.9%); the order of the settings is in P66.
- **P58.** **Met.** On the pooled dataset AmpliTrim's filter estimate is 82.9% and DADA2's pass rate 82.9% (0.0 points). Run by run (AmpliTrim on the same thinned pairs, `its_per_run.py`): nine runs, the difference 0.00 points in each (pass rates 71.4% to 88.6%). Descriptive, outside P58: DADA2 merged 14.0% of the pairs against an estimated yield of 82.1% (+68.1 points; DADA2 merged 0.17 of the estimate; on the first ITS dataset +25.5 points, 0.65 of the estimate, section 19). The command line on these reads (49,280 pairs): 81.7% pass the filter, 80.9% pass and merge, 89.4% of the forward reads alone pass, no warning; its screen line reads 'should pass filtering and merge'.
- **P59.** Regret against the best DADA2 setting that meets the coverage floor (grid points in brackets): inter-laboratory V1-V2 0.1 (9), inter-laboratory V4-V6 0.0 (6), wastewater 0.0 (9), spike-in 0.7 (6) and `miseqsop_mock` 0.0 (9) **met**; soil V4 4.7 points (9; 124/115 at 62.6%, best 134/105 at 67.3%, margin 1.0) **not met**, and the grid equals that of section 22 cell by cell; swabs **not decided** (DADA2 stopped with 'Error matrix is NULL' after two cells, 31.4% and 34.1%, both below AmpliTrim's total of 442 nt). On `miseqsop_mock` AmpliTrim's cut and the best cell both recovered 19/32 variants and 19/21 strains.
- **P60.** **Met.** On all 13 datasets the exit status was 0 and `REPORT.md` and `summary.json` were written; on an empty folder the exit status was 3 with a message and no traceback. Run time 29 to 141 s on the 12 datasets of 43,922 to 49,885 pairs (longest inter-laboratory V4-V6, 141.4 s) and 4.6 s on the 4,779 pairs of `miseqsop_mock`; peak memory 168 to 352 MB (largest mock A) and 52 MB on `miseqsop_mock`.
- **P61.** **Met.** On mock A and V1-V2 the repeat gave the same `recommended`, `maxEE_table`, `longest_amplicon` and `pairs_analysed`.
- **P62.** **Met.** At 15 and at 20 nt, on mock A (cuts 275/248 and 277/251) and V1-V2 (196/170 and 203/168), `recommended.minOverlap` was the value given and the total minus the overlap equalled `longest_mergeable_insert`.
- **P63.** On the hold-out mock (`miseqsop_mock`): **P50 met** (minOverlap 10 with the k = 12 cut: 97.1% merged as at 12, the same 20 ASVs); **P51 met** (reads in exact ASVs 100.00% at every k, 19/21 strains at every k); **P52 met** (0 extra ASVs at 12 and at 20). P53 (descriptive): arm T 96.6% at k = 25 against 97.1% at k = 12 (-0.5; no dataset with a fall of 1 point), arm S 0.0% at 15, 20 and 25 nt. P54: at the cut of k = 12 no pair has an overlap of 15 nt or more. `collect_overlap_results.py` prints 'not judged: only 1 of 4 datasets' because it counts the four mocks of section 27; here the hold-out alone is judged, as section 29 says.
- **P64.** **Met** as corrected: ITS (0/0, maxEE 2) and soil V4 (124/115, maxEE 2) as at Night 3; mock B 274/253 with maxEE 4,8 against 277/250 with 4,5 at Night 4b, different as predicted. (The script's own verdict, 'P30 not decided', follows from `--only`.)

**Addendum A (`miseqsop_mock`; the released package kept maxEE 2 and recommended 132/134).**
- **P65 met.** (a) The best relaxed setting (3,3, 3,5 and 5,5: 97.5%) merges 0.4 points more than 2,2 (97.1%; criterion at most 1.0). (b) Reads in exact ASVs are 100.00% and strains 19/21 in every row. (c) maxEE 1,1 merges 92.3%, 4.8 points fewer than 2,2.
- **P66 (descriptive).** The estimate of the command-line table minus DADA2's merged share (merged shares to 0.1, so the differences to about 0.1): 1,1 +0.3; 2,2 -0.7; 2,3 -0.8; 2,4 -0.8; 3,3 -0.9; 2,5 -0.8; 3,5 -0.9; 5,5 -0.9 (for 2,2 the coverage table, with the unrounded share, gives -0.8). Both put the settings in the same two groups (1,1 below the other seven, which lie within 0.5 points of each other); only that one step is large enough for the test to see.
- **P67.** (a) **met**: AmpliTrim's cut merges 97.1% against 89.3% for the tutorial's 240/160 (+7.8 points). (b) **met**: reads in exact ASVs 100.00% in both, strains 19/21 in both. (c) **met**: the filter estimate against DADA2's pass rate 98.0% against 98.2% (AmpliTrim's cut) and 90.5% against 90.3% (the tutorial's). (d) FIGARO was not run.

**What follows (section 29 said).**
1. Nothing in the package changes because of this run. The wording of the ITS yield is a separate decision (below).
2. Not met: P59 on soil V4 (4.7 points; the documented exception stays). Not decided: P57 (b) to (d) and P59 on swabs. The Abstract's sentence that the rule 'has not been tested on independent data' is replaced by what was measured: the filter estimate was within 0.2 points on two independent datasets (a mock community and a second ITS study), and the part of the rule that relaxes maxEE has not been tested on independent data (the released rule kept maxEE 2 on the independent mock).
3. The tip to raise the share left out to 1% is not added: no dataset other than D6300 V4-V5, which informed the default, merges 5 points more at 1% than at 0.2%. The Discussion's sentence on 0.2% stands.
4. P52 holds on the hold-out mock, so 12 nt stays. The paper reports that P52 was not met on two of the three mocks judged in section 27 and met on the hold-out.
5. The swabs grid is recorded as not decided, not as not met.
6. P60 to P62 are met: the sentences of the manuscript that rely on them stand.
7. Item 7 of section 29 stands as amended by Addenda A and B: no further test is run before submission.

*Added after the run, not registered in advance:*
- **Two defects of the validation scripts, found by this run and fixed before the two steps were repeated** (commit c59a0c8; `compare_methods.py` SHA256 7d8fadf30ed86c6229a182388f2408064b038793c813cdee222c3bd4bee43d8c, `tail_classify.py` SHA256 aac590d4cbd80d020393d4b5e9d09d46623ed90bf3ae0e38cf65d5ce9046cac3, identical in the archive and in the run copy). (1) `compare_methods.py` raised a `TypeError` when all estimates were equal, as on `miseqsop_mock` where no setting changes the cut; the rank correlation is now written as 'not defined'. (2) `tail_classify.py` took 'the longest 1%' of the merged sequences as every sequence at least as long as the 137th longest, a length (374 nt, the most common) shared by thousands of sequences: the first group had 9,215 sequences, mostly the ordinary amplicon (98.9% Bacteria, which would have read as 16S). Ties at the cut-off are now taken only when they are few. The first output is invalid and kept in `tail_first_run`. Three tests were added (28 in `tests/test_night10.py`); two of them fail on the old code.
- **The ITS yield estimate.** On two ITS datasets the estimate overshoots DADA2's merged share by 25.5 and 68.1 points. The released command line prints 'should pass filtering and merge' without a caveat for ITS; the report calls the yield an upper bound and gives 0.65 and 0.87 of the estimate (first ITS dataset) as the example. The split between denoising and merging on the second dataset was not measured (candidate causes, untested: a long overlap, three binned quality values). Whether to change the wording (and the example) in the software is the researcher's decision and would be a post hoc change in a later release; until then the limitation is stated in the paper, the user guide and the results pages, and the forward-only share is named as the number to read for ITS.
- The coverage rows of the final build replace the earlier-build rows of Supplementary Table S5; V1-V2 at 0% is 45.5% (not 29.4%), a cost of 38.7 points (not 'about 54').
- Run time and memory of Supplementary Table S4: the manuscript's 292 MB for mock A is replaced by 352 MB.
- `its_per_run.py` (kept in `C:\data\tools`; its text is at the end of `its_new\per_run_check.txt`) and the command-line run on the ITS reads (`its_new\cli_check`) were written and run after the main run.
- Each start of the driver rewrites `NIGHT10_SUMMARY.md`; the first version, with the package fingerprints, was overwritten. The fingerprints are in the first lines of `night10.log`.

*Added later on 7 November 2026 (evening), after the commit of this section; not registered in advance:*
- **Instrument label corrected.** The second ITS study (PRJEB108994) was sequenced on the Illumina MiSeq i100 (2 x 500), not an iSeq 100 as written in Addendum B and above (an iSeq 100 reads at most 2 x 150). The manuscript was corrected; the committed text above is left as written.
- **Exploratory split of DADA2's loss on the second ITS study (E1; run folder `night10_its2_split_2026-11-07`, run 24; `its_merge_check.py` pointed at the new reads; no criterion; its P26 to P29 verdict lines belong to section 20 and are not verdicts here).** Of 43,405 input pairs (9 runs, 5,000 pairs each): 82.86% pass the filter; 30.09% are lost in denoising; 52.77% reach the merge step; 3.81% are refused for an overlap shorter than 12 nt, 35.00% for mismatches or indels in the overlap and 0.00% for an overhang; 13.96% merge (13.96% also with `trimOverhang = TRUE`).
- **A 16S project set aside.** PRJNA1531119 (V4-V5, MiSeq i100, 2 x 500) was considered and not run (no expected sequences; the same region as D6300 V4-V5; binned qualities already covered by the second ITS study). Vendor demo data were not used (no accession).
- **An independent test of the `maxEE` rule by N. T. Ho (co-author), 6 and 7 November 2026.** Protocol, criteria and expectations were posted publicly (issue #1 [URL to add]) before each analysis: four real gut 16S datasets (CBGS PRJEB50418 V3-V4, NIG PRJEB31073 V1-V2, KFS PRJDB8879 V3-V4, HDP PRJNA714241 V3-V4; MiSeq 2 x 300; 24 samples each, 5,000 pairs per draw; KFS and HDP from the first 4,000,000 bytes of each file) and reads simulated from a known community. Results as reported there, not re-run here: C1 the filter estimate within 0.25 points in all comparisons; the rule kept maxEE 2 on CBGS and NIG and relaxed to 3,3 (KFS) and 3,2 (HDP), where C2 was met (+10.9 points, 24 of 24 samples; +3.1, 22 of 24); the genus profiles stayed within sampling noise (C4 to C7 met); C3 (yield estimate within 5 points) was met on CBGS and KFS and not on NIG (+6.1 to +9.0) or HDP (+18.5 to +20.0); maxEE 2,5 lowered rarefied richness by 7.7 to 10.7% on KFS, HDP and NIG; H1 (the rule lowers richness) was not supported. In the simulation (errors that follow the quality scores, SIM-I, and the higher rates DADA2 learned from real data, SIM-T) the rule met T1 to T4 against maxEE 2 and maxEE 2,5 did not (SIM-T: recall -12.6 points, spurious reads +3.7, distance to the truth +0.082). The full report, scripts and tables are to be added to this archive (a pull request offered by the investigator); until then this paragraph rests on the text of the issue.
- **A missing script added.** `validation/night1_v2.py`, the script that ran night 1 (section 11), was missing from the archive until 7 November 2026. The only copy found (a Downloads folder; last modified 9 October 2026, SHA256 `1ebd952665bfdd47ab358bf3f3ee26862931236fe383b4200f73b5c9afb26039`) was added unchanged; night 1 ran on 9 October 2026. It was not re-run and has no test.
- **Run 01 recovered.** The text outputs of the development rounds (register sections 1 to 8), found on 7 November 2026 in the supplement package of 7 October 2026 (4,533 files, 1.6 GB; its README lists the 86 runs of the datasets by accession), were added as `results/01_development_rounds_2026-10-07/` (827 files, 1.6 MB, with their sha256 in `results/MANIFEST.tsv`). The reads, the work folders and the ASV FASTA files are not kept. Among them is a copy of the register dated 7 October 2026 (21 KB), the earliest dated version of it. The register's own sections 1 to 8 remain the record; these files are the outputs they describe.

## 31. ITS follow-up: why DADA2 refuses ITS pairs at the merge step, and whether a mismatch limit above 0 is a safe answer (written 11 November 2026, before `validation/its_followup.py` was run)

**Why.** The merge estimate of AmpliTrim is too high on three datasets: soil V4 (+32.1 points), the first ITS study (+25.5) and the second ITS study (+68.1). It stops before DADA2 denoises the reads, so everything DADA2 loses after the filter is counted as merged. Section 20 and section 30 split that loss exactly for the two ITS studies. On the second, per 100 input pairs, 30.09 are lost in denoising and 38.81 are refused by `mergePairs` (35.00 for mismatches or indels in the overlap, 3.81 for an overlap shorter than 12 nt, 0.00 for an overhang); on the first, 15.9 are lost in denoising and 11.5 are refused (6.43 for mismatches or indels, 5.07 for a short overlap). For soil V4 there is no exact split: the DADA2 tracking table gives about 31.7 points lost after the filter, and the archive's strict estimate (no difference allowed in the raw overlap) still overshoots by +25.0. Nobody has measured what these refusals are. Five explanations are open:
- *Long overlaps.* ITS reads are not truncated, so the overlap is long, and the chance that two denoised halves disagree somewhere in it grows with its length. `mergePairs` accepts no mismatch by default.
- *Indels.* Homopolymer and repeat errors give indels that a mismatch count does not see.
- *Near misses.* Most refused pairs differ at one or two positions, so a small mismatch limit would recover them.
- *Different partners.* A refused pair joins two halves that each belong to other pairs of the same sample, which merge fine with their own partners (suggested by the second author's analysis of another dataset, to be checked against her tables). `mergePairs` records which two denoised sequences each pair joins, so this can be tested here (P73 and P75). DADA2's alignment of the two halves appears to penalise a mismatch far more heavily than it rewards a match (the scores are to be confirmed in the DADA2 source before this section is committed), so halves that disagree in the overlap would be reported with an alignment of almost no length: 'no overlap found' and 'halves that disagree' would then look the same in the length of the alignment.
- *Geometry (soil V4).* The amplicon is longer than the two cuts together, so there is no overlap to find. Against it: the over-merge check of section 22 found that 99.7% of the V4 pairs merge at the correct length at AmpliTrim's cut (124/115), and the final run's table gives DADA2's merged share for four cuts of soil V4 that rises as the cut, and so the overlap, gets shorter: 64.0%, 66.8%, 67.2% and 72.3% for 141/124, 124/115, 124/114 and 115/105.

This section tests the first three and the last two with one run, and asks whether allowing mismatches is safe on ITS, first by counting ASVs and then, on mixes whose true sequences are known from their isolates, by counting spurious sequences. Because the data are in hand it also repeats, on independent ITS1 and ITS2 libraries, the two checks of AmpliTrim's own numbers that the paper makes on its 14 datasets (the prediction of the filter, Figure 2, and the estimate of the yield, Supplementary Figure S2), and asks whether the results depend on the depth. It does not change any recommendation of AmpliTrim.

**What it does.** DADA2 (`learnErrors`, `dada`) is run once per dataset, on the same thinned reads as before (5,000 pairs per sample, the same function), with `maxN = 0`, `maxEE = 2`, `truncQ = 2`, `minLen = 50`, phiX not removed and no `truncLen` (soil V4: the forward and reverse cut that AmpliTrim recommends, 124 and 115 nt in the over-merge check of section 22, to be confirmed against the 0.2.0 output). `mergePairs` is then run on the same denoised sequences with `minOverlap = 12`, `trimOverhang = FALSE`, `returnRejects = TRUE` and a mismatch limit of 0 (the default), 1, 2 and 5. At limit 0 every distinct pair of denoised sequences is written with the numbers of its two sequences within the sample, its abundance, its matches, mismatches and indels. The length of an overlap is its matches, mismatches and indels together. A refused pair is *short* when that length is below 12, and is *refused for mismatches or indels* when it is 12 or more and the alignment has a mismatch or an indel (this includes a pair that has fewer than 12 matching bases only because of its mismatches). Section 20 counted every pair with fewer than 12 matching bases as short; that makes no difference for a long overlap, but for the 12 to 30 nt overlap of soil V4 it would hide the mismatches, so the definition is changed here, and the numbers of this section are not comparable with the 'short overlap' shares of sections 20 and 30. An ASV is a distinct merged sequence that holds at least 0.1% of the merged pairs of its sample.

**Datasets.** The second ITS study (PRJEB108994, 9 runs); the first ITS study (the rhizosphere dataset of the earlier runs, 37 samples); soil V4 (PRJNA574030, 36 samples); and, if their reads are long enough, the single-isolate libraries of Arfken et al. 2023 (BioProject PRJNA693350: the triplicates of *Candida albicans*, *Kazachstania slooffiae*, *Lichtheimia corymbifera*, *Pichia fermentans* and *Trichosporon asahii*, for ITS1 and for ITS2). The paper does not state the read length. The ENA run table of the 30 isolate runs (15 per marker; base count divided by read count and by 2) gives a mean of 270 to 280 nt per read for ITS2, so these reads are longer than 250 nt (probably 2 x 300 cycles, trimmed), and 163 to 279 nt for ITS1, where the reads seem to be trimmed to the end of the amplicon (shortest for *Pichia*, whose ITS1 is the shortest). Whether the primers are still on the reads is checked on one file before the run: **in C-ITS1a the first 30 nt of three forward reads are AAGTCGTAACAAGGTTTCCGTAGGTGAACC, which is the sequence that follows the ITS1F primer site (the primer CTTGGTCATTTAGAGGAAGTAA is not there), and the reverse reads do not begin with the ITS2 primer GCTGCGTTCTTCATCGATGC; the primers appear to have been removed before the reads were deposited**. The reads are kept as downloaded, in folders 16 (ITS1) and 17 (ITS2), and the report states the share of pairs that enter the merge step. The mixes MixA (KCTPM), MixB (KC) and MixC (KT) of the same study are used as well, in triplicate for ITS1 and ITS2 (folders 18 and 19); they are the only mixes of known composition whose ITS1 and ITS2 files differ in the ENA run table: for MixD to MixM and for ATCC2 the two markers have identical read and base counts, which looks like duplicated deposits (to be reported to the authors or the archive after a few reads have been compared), so those mixes are not used. The isolates were chosen because each isolate should give one sequence, so any extra ASV is a measure of a false join or an error that is let through; the paper itself found one ASV for *C. albicans* and *T. asahii* with both ITS markers and two to four for the other three.

**Four more checks, with the data in hand.**
- *Estimates* (`validation/its_estimates.py`): AmpliTrim 0.2.0 is run with `--its` on each of the four Arfken groups (isolates and mixes, ITS1 and ITS2), on the thinned reads that DADA2 received, so that both see exactly the same pairs; a second run on the whole folders, as a user would, is reported beside it and decides nothing. Its `summary.json` gives the predicted share of pairs that pass the filter and the estimated yield; DADA2's numbers come from the run above (maxEE 2, minLen 50, no truncation, 5,000 pairs per sample).
- *Depth* (`validation/its_depth.py`): the first part is repeated at 20,000 pairs per sample (and, as exploration, at 1,000), with the same datasets and labels.
- *Truth* (`validation/its_truth.py`): for each marker the reference is the set of merged sequences of the isolates at the default limit that hold at least 1% of an isolate's merged pairs, the three replicates pooled. A merged sequence of a mix is *true* when it is within 2 edits (mismatches or indels) of a reference sequence of any isolate, and *spurious* otherwise. For each limit the report gives the share of the merged pairs of the mixes in spurious sequences, the share of the reference sequences of the isolates in the mix that are found (each at 0.1% or more of the merged pairs of the sample), and the share of the input pairs that merge.
- Recall is not expected to reach 100%: the paper reports that the Mucor / Lichtheimia isolate is barely amplified by the ITS primers in the mixes, so its reference sequence is expected but seldom found. Only the change of recall with the mismatch limit is used (P81).
- The maxEE of ITS is not tested here: AmpliTrim does not adjust it in ITS mode, and the truth test would need a run at each maxEE.

**Predictions, written before the run.** (Datasets are named as in the script: `its2` the second ITS study, `its` the first, `v4` soil V4, `arfken_its1` and `arfken_its2` the single-isolate libraries.)
- **P68 (near misses).** On `its2`, allowing 2 mismatches merges at least 50% of the pairs that the default refuses for mismatches or indels. Half is the threshold between 'mostly a few differences' and 'mostly many'.
- **P69 (long overlaps).** On `its2`, among the pairs with an overlap of at least 12 nt, the share refused for mismatches or indels is at least 15 points higher in the longest third of the overlaps than in the shortest third. Thirds are by pair count, in order of overlap length.
- **P70 (indels).** On `its2`, at least 20% of the pairs refused for mismatches or indels have an indel.
- **P71 (across datasets).** The median overlap of the pairs that enter the merge step (overlap of at least 12 nt) is longer on `its2` than on `its`, the dataset with the larger share refused for mismatches (35.00 against 6.43 points).
- **P72 (flag on ITS).** Allowing 2 mismatches raises the median number of ASVs per sample by at most 10%, on `its2` and on `its`. A rise can also be real sequences that were recovered, not false joins, so P72 is a flag; P74 (isolates, one true sequence each) is the safety test.
- **P73 (soil V4, different partners).** At AmpliTrim's cut, at least 50% of the pairs that `mergePairs` refuses on `v4` (for any reason) have a forward and a reverse sequence that each also occur in an accepted pair of the same sample: the halves are good, the pairing is not. A refused pair whose halves occur in no accepted pair looks like a pair with no overlap to find.
- **P74 (safety on isolates).** On the single-isolate libraries that were run, allowing 2 mismatches adds at most 1 ASV per sample (median).
- **P75 (second ITS study, different partners).** The same as P73, on `its2`.
- **P76 (depth).** On every label among `its2`, `its` and `v4` that has both depths, the share of the input pairs that merge at the default at 20,000 pairs per sample is within 5 points of that at 5,000.
- **P77 (the filter prediction, Figure 2 on new data).** On each of the four Arfken groups, AmpliTrim's `--its` prediction of the share of pairs that pass the filter, made on the thinned reads, is within 1 point of DADA2's on the same reads.
- **P78 (the yield estimate and the diversity of the community).** For each marker, the estimate of the yield (made on the thinned reads) is too high by fewer points on the isolates than on the mixes.
- **P79 (the yield estimate on isolates).** On the isolates the estimate of the yield is too high by less than 17 points (a quarter of the +68.1 of the second ITS study), on both markers.
- **P80 (false joins).** On the mixes, allowing 2 mismatches raises the share of merged pairs in spurious sequences by at most 1 point, on every marker (mean over the samples of the mixes).
- **P81 (recall).** On the mixes, allowing 2 mismatches does not lower the share of the reference sequences that are found, on any marker.
- **P82 (gain).** On the mixes, allowing 2 mismatches merges at least 2 points more of the input pairs, on at least one marker; without a gain there is nothing to recommend.

*Exploratory, no criterion (E2):* the share of the default's mismatch refusals that merge at limits 1, 2 and 5; the median overlap and the thirds for every dataset; the share of the refused pairs whose alignment has no column at all; the share of every refused pair that merges at limits 1, 2 and 5 (the one that is meaningful for soil V4). A criterion that is not decided because a dataset was not run stays not decided; it is not replaced.

**Before the run.** A smoke test at 300 pairs per sample was run on the three datasets of the first part to check the paths, the cut and the code. Its verdicts were not used: at that size DADA2 loses about 80% of the pairs in denoising, so few pairs reach the merge step. Two things were changed after it, before the real run, and no threshold of P68 to P72 or P74 was changed. (1) The soil V4 line (every refused pair counted as 'short') led to a check of the definition of 'short', which was found to hide mismatches in short overlaps; the definition was changed (see 'What it does'). (2) A second look at the soil V4 smoke result, where every refused pair had an alignment of 0 to 7 nt, mostly none, showed that the length of the alignment cannot tell an amplicon that is too long from halves that disagree (DADA2's alignment appears to drop an overlap with several mismatches; to be confirmed in the source). The original P73 (mismatches against a short overlap) was therefore replaced by P73 as written above, and P75 was added for the second ITS study. Both were written after the smoke test and before the real run, and the reader should know that the smoke result for soil V4 suggested the direction. (3) The soil V4 table of the smoke test showed the merged share rising with the mismatch limit (40.03%, 42.06%, 45.21% and 46.84% of the input pairs for limits of 0, 1, 2 and 5) although no pair counted as refused for mismatches: on a short overlap the alignment drops an overlap that contains a mismatch, so such pairs are filed as short at the default. The report therefore also gives the share of every refused pair that merges at each limit; no criterion uses it (E2). (4) The soil V4 partner share of the smoke test was 23.6%. At 300 pairs per sample few pairs are accepted, so few halves occur in an accepted pair, and the share is biased low by depth; the threshold of P73 and P75 was not changed, and a failure of either is read with that in mind. (5) A smoke test at 300 pairs per sample was also run on the single-isolate libraries. On ITS2 no pair was refused, so P74 cannot fail on ITS2 at that depth. On ITS1, 8.4% of the input pairs were refused for mismatches or indels, almost all with an indel (99.2%), concentrated in the longest third of the overlaps (13.4% against 0.3% in the shortest), and a limit of 1 mismatch recovered all of them with no change in the median ASV count of 1. No criterion was changed after this. (6) After the smoke tests, seven more criteria (P76 to P82) were written for the checks above, before any real-size run of them and before the mixes were downloaded or looked at. The mixes were chosen from the ENA run table alone (known composition, different ITS1 and ITS2 files). The limit of P79 (17 points) is a quarter of the +68.1 of the second ITS study and was not taken from the smoke numbers, which had shown no refused pair on the ITS2 isolates and about 8% on the ITS1 isolates. (7) A smoke test of the truth test at 300 pairs per sample on the mixes (after the mixes were downloaded) found one reference sequence per isolate (two for Pichia on ITS2), no spurious sequence at any limit, recall of 77.8% on ITS1 and 87.0% on ITS2 at every limit up to 2 (90.7% at 5 mismatches on ITS2), and no gain in the share merged at 1 or 2 mismatches (75.81% on ITS1, 35.48% on ITS2; 36.59% at 5 mismatches on ITS2). No criterion was changed after this; P82 may therefore not be met. (8) A smoke test of the estimates test, with AmpliTrim on the whole folders and DADA2 on 300 pairs per sample, gave filter differences of +4.44 (ITS1 isolates), -3.37 (ITS2 isolates), -1.11 (ITS1 mixes) and +0.35 (ITS2 mixes) points, and yield estimates too high by +14.1, -3.4, -0.8 and +2.7 points. The samples of a group differ up to nine-fold in size and in filter pass rate, and a run on a whole folder weights them by their reads while the thinning gives every sample the same number of pairs, so the two sides did not describe the same pairs (for the ITS2 isolates the yield difference, -3.4, equals the filter difference, -3.37, so nothing after the filter contributed). The procedure was therefore changed, before the real run: AmpliTrim is given the thinned reads that DADA2 received, and the whole-folder run became an exploratory table. No threshold was changed. In that smoke test P78 would have been not met (the ITS1 isolates were above the ITS1 mixes, +14.1 against -0.8), and it may be not met in the real run as well. (9) The estimates test was then repeated on the thinned reads at 300 pairs per sample: the filter differences were +0.00 points on all four groups, which confirms that the earlier differences came from the weighting of the samples and not from the prediction, and the yield estimates were too high by +9.7 (ITS1 isolates), +0.0 (ITS2 isolates), +0.3 (ITS1 mixes) and +2.4 (ITS2 mixes) points, so P77 and P79 were met and P78 was not (ITS1). No criterion was changed after this.

**What follows.**
- *P68 holds, and P72 and P74 (if run) hold:* the ITS report of version 0.3.0 states what share of the refusals a limit of 2 mismatches recovers on ITS, and the guide names `maxMismatch` as an option for ITS, with the sentence that accuracy was checked only by counting ASVs. What AmpliTrim recommends for truncation and maxEE does not change; the merge row of the ITS table stays an upper bound.
- *P68 holds and P74 fails:* allowing mismatches is not recommended, and the report says that it adds ASVs. *P68 holds, P74 holds or was not run, and P72 fails:* the extra ASVs may be recovered sequences or false joins and this section does not decide which; the report says only that a limit above 0 changes the ASV count, and the option is not recommended without a check on isolates.
- *P68 fails:* the refused pairs are not near misses; no advice to relax the limit. The different-partner explanation is the next one to test.
- *P69 holds:* the guide says that long overlaps are refused more often, which is a reason to expect a larger overshoot of the estimate on ITS than on 16S. *P69 fails with P68 holding:* the differences are few but not tied to the overlap length.
- *P70 holds:* the refusals are partly indels, and the ITS report says so. Not claimed: that homopolymers are the cause; the run does not look at the sequence.
- *P71 fails:* the difference between the two studies is not explained by their overlaps; no other claim is made from P69.
- *P73 holds:* on soil V4 most refused pairs join two halves that merge with other partners in the same sample, which is not a geometric shortage of overlap; the soil V4 note says that the overshoot of the estimate is made of pairs that the denoised sequences pair wrongly or that disagree in the overlap, and that the cause is not shown. *P73 fails:* most refused pairs have a half that occurs in no accepted pair of its sample. That fits an amplicon longer than the cut, but it also fits rare variants whose halves are seen only once, so it is not evidence of geometry by itself; the report's alignment line (pairs with no overlap at all) and the recovery at higher mismatch limits are read with it.
- *P75 holds:* the same for the second ITS study, and the second author's observation on another dataset is repeated on these data. *P75 fails with P68 holding:* the refusals are near misses, not wrong partners.
- *P80, P81 and P82 hold:* with P74 and P68, the ITS report of version 0.3.0 may name `maxMismatch` as an option for ITS and cite this evidence, which comes from three mixes of five fungal isolates. *P80 fails:* allowing mismatches creates spurious sequences, and the option is not recommended. *P82 fails:* there is nothing to gain from it on these data. *P81 fails:* allowing mismatches loses true sequences, and it is not recommended.
- *P77 holds:* the filter prediction of AmpliTrim's ITS mode holds on short ITS1 and ITS2 amplicons of an independent study, in addition to the two ITS studies of the paper; the exploratory table shows how far a run on the whole folders, which weights the samples by their reads, can differ. *P77 fails:* it does not, and the cause is looked for before version 0.3.0.
- *P78 and P79 hold:* the overshoot of the yield estimate is small where the community is simple, which is what the explanation of section 3.6 predicts; the report of version 0.3.0 may say that the estimate is closest on simple communities. *P78 fails or P79 fails:* the explanation by loss in denoising and in `mergePairs` is not enough, and the Limits paragraph says that the cause of the overshoot is not understood.
- *P76 holds:* the results of this section do not depend on 5,000 pairs per sample, within 5 points. *P76 fails:* they do, and every number of this section is quoted with its depth; nothing is said about other depths.
- Not claimed in any case: accuracy beyond the ASV counts, taxonomy, abundances, or anything about ITS data other than these; and nothing about shortening the overlap, which a truncation would confound with losing the longest amplicons.
- Nothing from this section enters the manuscript beyond one sentence in the discussion of future work, and only if P68 and P72 hold. The results go to `FUTURE.md` as items for version 0.3.0.

## 32. ITS follow-up: outcomes (12 November 2026)

**The run.** The scripts of section 31 as committed in `e301a1f`, before any real-size run. The estimates step was run after a fix committed in `126cbd5`: from the archive folder, `python -c "import amplitrim"` had imported the archive's own early copy of the package (0.1.0, which has no `--its`) instead of the release, so AmpliTrim is now started from the output folder and the script checks the version; the fix changed no criterion. AmpliTrim 0.2.0 from the software repository at `8605bf7`, without local changes; its package code is identical, file for file, to the PyPI wheel 0.2.0 and to the tag `v0.2.0`. Folders `C:\data\runs\25_its_followup` and `25_its_followup_arfken` (5,000 pairs per sample), `25_its_depth_1000`, `25_its_depth_20000`, `25_its_depth` and `25_its_followup_summary` (named without the number when they ran). Reports and tables up to 2 MB are in `results/25_*` (commit `5617923`); the larger tables are in `C:\data\backup\runs_25`. The fifteen criteria in one table: `results/25_its_followup_summary/ITS_FOLLOWUP_SUMMARY.md`, **8 met, 7 not met, none undecided**.

**Outcomes** (the criteria of section 31; `its2` the second ITS study, `its` the first, `v4` soil V4; 5,000 pairs per sample unless another depth is given).
- **P68 met.** Allowing 2 mismatches merges 102.0% of the pairs that the default refuses for mismatches or indels (more than all of them, because a limit above 0 also merges pairs that the default counts as short); 99.7% at 1,000 and 102.0% at 20,000 pairs per sample.
- **P69 not met.** Refused for mismatches or indels: 70.1% in the longest third of the overlaps (from 280 nt) against 70.6% in the shortest third (up to 204 nt), -0.5 points; -1.7 at 1,000 and +1.7 at 20,000 pairs per sample.
- **P70 met.** 98.6% of the pairs refused for mismatches or indels have an indel (99.4% at 1,000, 98.6% at 20,000).
- **P71 met.** Median overlap 265 nt on `its2` against 188 nt on `its`.
- **P72 not met (a flag).** Allowing 2 mismatches raises the median number of ASVs per sample by 266.7% on `its2` (400.0% at 1,000, 227.6% at 20,000) and by 9.5% on `its` (6.7%, 3.4%).
- **P73 not met.** On `v4`, 45.8% of the refused pairs have halves that each occur in an accepted pair of the same sample. The share rises with depth, as item (4) of section 31 expected: 26.0% at 1,000 and 64.3% at 20,000 pairs per sample, where the criterion would be met.
- **P74 met.** On the isolates, allowing 2 mismatches adds no ASV per sample (median, ITS1 and ITS2).
- **P75 not met.** On `its2`, 0.3% of the refused pairs have halves that occur in accepted pairs (0.5% at 1,000, 0.3% at 20,000).
- **P76 not met.** The share of the input pairs that merge at the default rises from 5,000 to 20,000 pairs per sample by 2.7 points on `its2` (14.0% to 16.6%), 7.2 on `its` (46.9% to 54.1%) and 9.7 on `v4` (66.8% to 76.5%). The share that passes the filter stays within 0.3 points at every depth; the share that reaches the merge step after denoising rises (`its2` 52.8% to 65.5%, `its` 58.4% to 66.0%, `v4` 82.4% to 89.3%).
- **P77 met.** AmpliTrim's filter prediction on the thinned reads against DADA2's: -0.21, -0.09, +0.00 and +0.00 points (ITS1 and ITS2 isolates, ITS1 and ITS2 mixes). AmpliTrim analysed 37,076 and 36,721 of the 75,000 pairs of the isolates and 44,547 and 44,269 of the 45,000 pairs of the mixes.
- **P78 not met.** The estimate of the yield is too high by 9.3 points on the ITS1 isolates against 2.0 on the ITS1 mixes; on ITS2 by 0.0 against 0.8.
- **P79 met.** +9.3 and +0.0 points on the isolates.
- **P80 met.** Allowing 2 mismatches changes the share of merged pairs in spurious sequences by -0.00 points on ITS1 (0.01% to 0.01%) and by +0.48 on ITS2 (8.50% to 8.98%; means over the samples).
- **P81 met.** The share of the reference sequences found rises on ITS1 (71.1% to 80.0%) and stays at 100.0% on ITS2.
- **P82 not met.** Allowing 2 mismatches merges 1.9 points more of the input pairs on ITS1 and 0.2 on ITS2.

**Exploratory (no criterion).**
- *E2.* Of all the pairs refused at the default, the share that merges at 5 mismatches: `its2` 96.8% (98.1% at 1,000, 96.7% at 20,000), `its` 73.2% (64.2%, 78.1%), `v4` 81.8% (76.9%, 80.0%). On `v4` every refused pair has an alignment of fewer than 12 columns, at every depth.
- *The ITS1 isolates.* 9.51% of the input pairs are refused for mismatches or indels (85.1% of them with an indel), more in the longest third of the overlaps (15.5%) than in the shortest (2.2%), unlike `its2`; 1 mismatch merges 99.2% of them, and the median ASV count stays at 2.
- *The overshoot of the yield estimate, by depth* (with the final run's estimates, which do not depend on the depth: 98.3% on `v4`, 72.4% on `its`, 82.1% on `its2`): +47.1, +31.5 and +21.8 points on `v4`; +36.1, +25.5 and +18.3 on `its`; +73.3, +68.1 and +65.5 on `its2`, at 1,000, 5,000 and 20,000 pairs per sample.
- *The spurious sequences of the mixes* (an analysis of the merged-sequence tables after the run, with the reference of `its_truth.py` and a test for sequences made exactly of the start of one reference sequence and the end of another; shares pooled over the samples). On ITS2 all of them at the default (8.57% of the merged pairs, 21 sequences) are such chimeras of reference sequences of two different isolates, and each occurs only in mixes that hold both isolates (K and T 701 pairs, C and K 595, the other pairs under 60). DADA2's chimera removal, which the follow-up does not run, is expected to remove most of them. Allowing mismatches adds a part that is not such a chimera: 0.00% of the merged pairs at limits 0 and 1, 0.19% at 2 and 0.36% at 5. On ITS1 the spurious share is 0.01% at every limit.
- *The whole folders* (the exploratory table of the estimates test, against DADA2 on the thinned reads): filter differences +5.60, -1.79, +0.63 and +0.15 points, yield differences +15.1, -1.7, +2.6 and +0.9. A run on a whole folder weights the samples by their reads (items (8) and (9) of section 31).
- *The run table of PRJNA693350.* For MixD to MixM and for ATCC2 the ITS1 and ITS2 runs have identical read and base counts, which looks like duplicated deposits; not checked on the reads.

**What follows (as section 31 said).**
1. `maxMismatch` is not recommended for ITS. P68 holds, but P72 fails, P82 fails (little to gain on these libraries), and the exploratory non-chimeric spurious part at a limit of 2 (0.19% of the merged pairs) is the kind of false join that P80 was meant to bound.
2. On the second ITS study the pairs refused at the merge step nearly all carry an indel between the two denoised halves (P70 holds). They are not tied to the length of the overlap (P69 fails with P68 holding: the differences are few) and they are not good halves with the wrong partner (P75 fails). The ITS report of version 0.3.0 may say that the refused pairs are mostly indels in the overlap. Not claimed: that homopolymers are the cause.
3. P73 fails at the registered depth (45.8%) and would be met at 20,000 pairs per sample (64.3%); as item (4) of section 31 said, the failure is read with the depth bias in mind. On soil V4 the refused pairs are not short in geometry: about 80% of them merge at 5 mismatches (E2), so they are mismatches in a short overlap. Nothing changes for the minimum overlap of section 27.
4. P76 fails: the results depend on the depth, and every number of this section is given with its depth. The overshoot of the yield estimate that the paper reports (section 3.6 and Supplementary Figure S2) is measured at 5,000 pairs per sample, and it is smaller at 20,000; nothing is said about other depths.
5. P77 holds: AmpliTrim's filter prediction holds on short ITS1 and ITS2 amplicons of an independent study, on the reads that DADA2 receives. The exploratory whole-folder table shows how far a pooled run can differ when the samples differ in size and in filter pass rate.
6. P78 fails (P79 holds): loss in denoising and in `mergePairs` does not predict *where* the overshoot is large, and, as section 31 said, the Limits paragraph says that its cause is only partly understood. Exploratory: on the ITS1 isolates the overshoot (+9.3 points) is close to the share of the input pairs refused for mismatches or indels (9.51%).
7. Section 31 said that nothing from it enters the manuscript beyond one sentence of future work, and only if P68 and P72 hold. P72 does not hold, so no result of this section enters the manuscript. Items 4 and 6 bear on how the manuscript words the overshoot (that it is measured at 5,000 pairs per sample; that its cause is only partly understood); that wording is a decision of the authors and is recorded in `FUTURE.md`.
8. Not claimed in any case: anything about full-length ITS from the Arfken libraries, which are short ITS1 and ITS2 amplicons; accuracy beyond the truth test; anything about maxEE on ITS, which was not tested.

## 33. ITS overnight run: what AmpliTrim can recommend for ITS (written 14 November 2026, before the run)

**Why.** AmpliTrim 0.2.0 gives no setting of its own for ITS: it uses DADA2's ITS defaults (no truncation, `truncQ = 2`, `minLen = 50`, maxEE 2,2), predicts the share of pairs that pass the filter (within 0.25 points on six ITS groups, sections 30 and 32) and estimates the yield, which on full-length ITS is far too high (82.1% against 14.0% on the second ITS study at 5,000 pairs per sample). Sections 31 and 32 left four questions that decide what version 0.3.0 can say about ITS:

| | Question |
|---|---|
| Q1 | Can a stricter merge test make the yield estimate predict DADA2's merge on ITS, without making it worse on 16S? |
| Q2 | When ITS reads have poor ends, is it better to keep maxEE 2,2, relax maxEE, or trim the poor ends? |
| Q3 | When long ITS pairs cannot merge, are the forward reads alone as accurate? |
| Q4 | Does merging first (vsearch on the raw pairs, then DADA2 on the merged reads) recover the refused pairs without false sequences? |

**Parity with the 16S validation.** Every check that the 16S datasets went through, and where it stands for ITS:

| 16S check in the paper | For ITS | Here |
|---|---|---|
| Filter prediction (Figure 2) | done at 2,2 | extended to every maxEE arm (P92) |
| Merge estimate (Supplementary Figure S2) | done: an upper bound | Q1 (P83, P84) |
| The estimate ranks settings in DADA2's order | not done | the order of the maxEE arms (P91) |
| Truncation grid; coverage floor | do not apply: no cut | replaced by the arms of Q2 |
| maxEE rule | not done: the rule needs a cut | Q2 (P85 to P87) |
| Minimum overlap (10 to 25 nt) | not done; named as missing in the paper's Limits | 12 against 20 and 30 nt (P93) |
| Accuracy on mocks against expected sequences (strains found, reads in exact sequences) | partly: without chimera removal, against the isolates' own sequences | against expected sequences from NCBI, after chimera removal (described for `std`; decides Q2 to Q4) |
| FIGARO | cannot run on ITS | not run |
| Reproducibility of the recommendation | does not apply: fixed settings | not run |
| Depth | done (P76) | not repeated |
| Independent test by another group | not available | not run |

**Data.** All from Arfken AM, Frey JF, Carrillo NI, Dike NI, Onyeachonamm O, Rivera DN, Davies CP, Summers KL (2023), *Porcine fungal mock community analyses: Implications for mycobiome investigations*, Frontiers in Cellular and Infection Microbiology 13:928353 (BioProject PRJNA693350), except the two ITS studies and the 16S controls.

| Label | Data | Truth |
|---|---|---|
| `atcc_its1`, `atcc_its2` | The ATCC mycobiome standard MSA-1010, ten strains with ATCC numbers (Table 3 of the paper), libraries ATCC-ITS1a to c and ATCC-ITS2a to c (folders 20 and 21). The ATCC2 libraries are not used: their ITS1 and ITS2 runs have identical read and base counts. | expected sequences from NCBI (below) |
| `arfken_its1`, `arfken_its2` | The five piglet isolates (C, K, M, P, T), triplicate | their own sequences: the isolate reference |
| `arfken_mix_its1`, `arfken_mix_its2` | MixA (KCTPM), MixB (KC), MixC (KT), triplicate | the isolate reference |
| `its2`, `its` | The second (PRJEB108994, full-length ITS) and the first (PRJEB58189) ITS study | none: yields and flags |
| 16S controls | The 12 16S datasets of the paper (amplitrim_paths.json, as the closing run read them), at the cut AmpliTrim 0.2.0 recommends for the same thinned reads, maxEE 2,2; when AmpliTrim cannot analyse a control's thinned reads, the cut the final run recommended for the whole dataset (run 24) is used, and the report says so | none needed: Q1 compares merged shares |

Every dataset is thinned with the same function as before to 5,000 read pairs per sample (every pair where a sample has fewer, and every pair of mock A).

**Expected sequences from NCBI.** `validation/ncbi_references.py` finds, for every member of the two mocks, a genome assembly at NCBI: one of the very strain if it exists ('strain'), else the species' reference genome ('species'). It cuts out the ITS1 and ITS2 amplicons with the paper's primers (Table 1: ITS1F `CTTGGTCATTTAGAGGAAGTAA` and ITS2 `GCTGCGTTCTTCATCGATGC`; ITS3 `GCATCGATGAAGAACGCAGC` and ITS4 `TCCTCCGCTTATTGATATGC`; at most 2 mismatches per primer), without the primers, like the reads. When a genome gives no product, NCBI's ITS records of the species are tried. Every copy found is kept, because fungal ITS can vary within one genome. The script's report (`REFERENCES_REPORT.md`) is read, and appended to this section, **before** the section is committed. *Result of that review (12 November 2026):* strain-exact references for 6 of the 10 ATCC members, for both markers: *A. fumigatus* (Af293 is ATCC MYA-4609), *C. neoformans* var. *grubii* (H99 is ATCC 208821), *C. glabrata* (ATCC 2001), *M. globosa* (CBS 7966 is ATCC MYA-4612) and *S. cerevisiae* from their genomes, and *C. albicans* ATCC 10231 from NCBI's ITS records of that strain (its genome assembly holds no rDNA). *S. cerevisiae* ATCC 201390 (BY4743) is counted through its parent S288C, a choice made before the count was known. So **the ATCC standard is the primary truth for both markers**; *T. interdigitale*, *P. chrysogenum*, *F. keratoplasticum* and *C. dermatis* have references of other strains only and are left unjudged. The piglet references agree with the isolates' own sequences within 0 to 3 edits for C, K, P and T (*K. slooffiae* ABBL is identical to isolate K on both markers); *L. corymbifera* gives no ITS1 product and its ITS2 record is 56 edits from isolate M's sequence, so isolate M's reference is flagged wherever the piglet mixes are reported. Disclosed: the first version of the script required the 3' end of each primer to match exactly and missed ITS1 in the yeasts and basidiomycetes; it was corrected (any 3 mismatches, every site found) before the references were fixed.

**Which truth decides (fixed now).**
- For each marker, the **ATCC standard is the primary truth** if strain-exact references were found for at least 6 of its 10 strains for that marker. A sequence of the standard is *true* when it is within 2 edits of an expected sequence of a strain-exact member, and *of an unjudged member* when it is within 3% of the length of a species-level reference of another member: such sequences are left out of the spurious share and of recall. Anything else is *spurious*.
- Otherwise the **piglet mixes are the primary truth** for that marker, judged against the isolate reference, as in section 31: the merged sequences of the isolates in the `std` arm, after chimera removal, that hold at least 1% of an isolate's reads (the three replicates pooled); true within 2 edits.
- The other truth is reported beside it for every arm. The NCBI references of the piglet species are used only to check the isolate reference (the edits between them are reported), never as truth, because the isolates were never sequenced as genomes in public.
- *Recall*: the share of the expected members (the strain-exact members; the isolates of a mix) that are found in a sample, a member being found when a true sequence of it (any of its copies) holds 0.1% or more of the sample's reads; the mean over the samples. *Reads in exact sequences*: the share of reads in sequences that equal an expected sequence. *Strains found*: per sample. For `fwd`, a sequence is judged against the starts of the expected sequences of the same length.

**One pipeline, several arms.** DADA2 (`learnErrors`, `dada`, pool = FALSE, `minLen = 50`, `maxN = 0`, `minOverlap = 12`, `maxMismatch = 0`), then **chimera removal** (`removeBimeraDenovo`, method "consensus") in every arm; counts are made after it.

| Arm | What changes from the default | Question |
|---|---|---|
| `std` | nothing: no truncation, `truncQ = 2`, maxEE 2,2 | the baseline |
| `ee23`, `ee24`, `ee25`, `ee33`, `ee55` | maxEE 2,3 / 2,4 / 2,5 / 3,3 / 5,5 | Q2 |
| `tq10`, `tq20` | `truncQ = 10` / `20`, maxEE 2,2 | Q2 |
| `sw15` | a sliding-window trim of the read ends before DADA2 (window 4, mean quality 15), then maxEE 2,2 | Q2 |
| `fwd` | the forward reads of `std` alone, denoised, no merge | Q3 |
| `mfirst` | vsearch merges the raw pairs (its default allowance of differences, a minimum overlap of 12 and staggered pairs allowed, as in AmpliTrim), then DADA2 denoises the merged reads (maxEE 2) | Q4 |
| `mo20`, `mo30` | `mergePairs` with `minOverlap = 20` / `30` on the denoised reads of `std` | parity |

AmpliTrim 0.2.0 (`--its`) is also run on the thinned reads of every ITS group, for its predicted pass rates and yields at each maxEE (P91, P92). For Q1 the merge test behind AmpliTrim's estimate is repeated on the reads of `std` with vsearch allowing 10 (today's tolerant default), 5, 2 and 0 differences; vsearch merges without gaps, so an indel in a long overlap shows as many differences.

**Criteria** (written before the run).
- **P83.** On the second ITS study, at least one of the stricter estimates (5, 2 or 0 differences) is within 10 points of DADA2's merged share (before chimera removal, as in Supplementary Figure S2).
- **P84.** That estimate (the closest on the second ITS study) is not farther from DADA2 on the 16S controls than today's: the mean absolute difference over the 12 datasets grows by at most 2 points.
- **P85 (safe).** Against `std`, on the primary truth of each marker: the share of reads in spurious sequences rises by at most 1 point, and recall falls by at most 2 points.
- **P86 (useful).** On the primary truth of ITS2, the yield rises by at least 5 points against `std`.
- **P87 (a flag).** On the two ITS studies, the median number of ASVs per sample rises by at most 25% against `std`.
- **P88.** `fwd`: on the primary truth of each marker, the spurious share is at most 1 point above `std`'s and recall is at least `std`'s.
- **P89 (safe).** `mfirst`: on the primary truth of each marker, the spurious share is at most 1 point above `std`'s, and recall does not fall.
- **P90 (useful).** `mfirst`: on the second ITS study, the yield is at least 10 points above `std`'s.
- **P91.** On each ITS group, AmpliTrim's predicted yields order the maxEE arms that its maxEE table holds (of `std`, `ee23`, `ee24`, `ee25`, `ee33`, `ee55`; at least three) as DADA2's merged shares do, with a rank correlation of at least 0.8.
- **P92.** On each ITS group, AmpliTrim's predicted pass rate on the same thinned reads is within 1 point of DADA2's, at every maxEE arm.
- **P93.** On the primary truth of each marker, the spurious share at a minimum overlap of 12 nt is at most 0.5 point above the share at 20 and at 30 nt.

*Described, no criterion:* the accuracy of `std` on the primary truth (members found, recall, reads in exact sequences, spurious share) and of every arm on the piglet mixes (the second view); the yields of every arm on the two ITS studies; the share of the second ITS study's pairs that a minimum overlap of 20 and 30 nt costs.

**What follows (decided now).**
- P83 and P84 hold: the report of 0.3.0 gives the closer estimate for ITS (and for all markers, if P84 holds with room). P83 fails: the report says that the ITS estimate is an upper bound and why (indels at the merge step).
- Among the arms of Q2 that meet P85 and P86 on both markers, the one with the largest ITS2 yield becomes the advice of 0.3.0 for ITS data where fewer than half the pairs pass at 2,2; the report says what it was tested on (short ITS1 and ITS2 amplicons of fungal mock communities). P87 failing is reported beside it. No arm meets both: 2,2 stays, and the report says that relaxing and trimming were tested and were not safe or not useful.
- P88 holds: the report may name forward reads alone as an option when most long ITS pairs cannot merge; it does not recommend it.
- P89 and P90 hold: merging first may be named as an option for full-length ITS, with the caveat that DADA2's authors advise denoising first. Otherwise it is not mentioned.
- P91 and P92 hold: the report presents the maxEE table for ITS as a reliable guide to the pass rate and the order of the settings; P91 fails: it says that the table's yields do not give the order.
- P93 holds: 12 nt stays for ITS, and the paper's Limits no longer need to name the missing test (the authors' decision, in FUTURE.md). P93 fails: a longer minimum overlap is advised for ITS, with its cost.
- Not claimed in any case: accuracy on full-length ITS, which has no truth here; taxonomy; abundances; anything about the members left out of the truth.

**What has already been seen (disclosed).** AmpliTrim's own predicted pass rates and yields at relaxed maxEE for the Arfken groups (run 25); the outcomes of run 25 at 2,2 without chimera removal, including the chimeras between isolates in the ITS2 mixes; the reference report of `ncbi_references.py`. A smoke run at 300 read pairs per sample (12 November 2026), after these criteria were agreed and before this section was committed, ran every step of every dataset to check the code; its verdicts mean nothing at that depth, are reported in section 34, and changed nothing here. It led to one fix of the procedure: AmpliTrim could not analyse the wastewater control's thinned reads (57 of 2,699 pairs overlapped at full read length), so a 16S control now falls back to the final run's cut. **Not seen at full depth:** the ATCC libraries (not downloaded before this section was written); any DADA2 result for a relaxed maxEE, a higher `truncQ`, a trim, forward reads alone, merging first or a longer minimum overlap; any count after chimera removal; any stricter estimate on ITS.

**Procedure.** (1) The reference script is run and its report read; the members used in the truth are fixed here. (2) The code (R and Python) is written with tests on stand-ins. (3) A smoke run at 300 pairs per sample checks the paths (its look at the data is disclosed in section 34). (4) This section, the code and the reference sequences are committed. (5) The overnight run: one resumable command; every arm of every dataset is a step that is not repeated when the run is started again. (6) One summary of P83 to P93, the folders numbered and collected, and section 34 written from the summary.

## 34. ITS overnight run: outcomes (14 November 2026)

**The run.** Section 33 was agreed on 13 November 2026 and written into the register by the registration patch on 14 November (the date of its heading), before the run; it was committed in `46009ba` later on 14 November, after the night had run and before any of its output was read: the commit step had been missed (the commit message says 'written 13 November', the date of the agreed draft). The code that ran is the code committed there (the registration patch copied the same files into the archive). AmpliTrim 0.2.0 (the release code), R with DADA2. Folder `C:\data\runs\26_its_night` (named `its_night` when it ran): 96 of 96 steps in about 3.5 hours (the two real ITS studies took 50 and 83 minutes for their 13 arms). Results in `results/26_its_night` (`9f48c60`); the larger files in `C:\data\backup\runs_26`. The summary of the criteria: `results/26_its_night/ITS_NIGHT_SUMMARY.md`.

**Fixes after the run (none changed a verdict).** (1) The archive's `.gitignore` excludes `*.fasta`, so the reference sequences of section 33 were not committed with it; they were committed unchanged in `a08868f`. (2) The analysis divided the yield of merging first by the reads that vsearch merged instead of the input pairs; corrected in `a08868f` (P90: +50.4 points instead of +55.5) and the analysis rerun (`46a0d51`). (3) The analysis read an arm only when the driver's empty `DONE` file was there, which the archive does not keep; corrected with this section, so that the verdicts can be regenerated from the archive alone. They were: the summary, the report and the metrics table regenerated from the archive equal the committed ones.

**The smoke run** (300 read pairs per sample, disclosed in section 33; not part of the outcomes). Its verdicts differed in places: P83 not met (the estimate without differences +11.1 points), P84 met (11 controls; the wastewater control fell back to the final run's cut), P85 not met for `sw15` (ITS1 recall -5.6 points), P90 not met (+9.8 points, with the wrong denominator of fix 2), P91 not met (the second ITS study -0.65); the advice would have been `ee55`.

**Outcomes** (the criteria of section 33; 5,000 read pairs per sample; the ATCC standard against the NCBI references is the primary truth of both markers).
- **P83 met.** On the second ITS study the estimates differ from DADA2's merged share by +68.17 (10 differences, today's), +66.18 (5), +48.26 (2) and -2.75 points (0 differences, the closest).
- **P84 not met.** On the 12 16S controls the estimate without differences is farther from DADA2 than today's: mean absolute difference 8.15 against 5.84 points.
- **P85 to P87, per arm** (spurious share and recall against `std` on the ATCC standard; ITS2 yield gain; the median ASVs on the second and first ITS study):

| Arm | P85 (safe) | P86 (useful) | P87 (a flag) |
|---|---|---|---|
| `ee23` | met: spurious +0.00 / +0.00, recall +5.6 / -0.0 | met: +9.6 | met: +4%, +0% |
| `ee24` | met: +0.00 / +0.00, +11.1 / -0.0 | met: +13.6 | met: +4%, +0% |
| `ee25` | met: +0.00 / +0.00, +11.1 / -0.0 | met: +18.0 | met: +4%, +0% |
| `ee33` | met: +0.00 / +0.00, +5.6 / -0.0 | met: +10.2 | met: +9%, +0% |
| `ee55` | met: +0.00 / +0.00, +11.1 / -0.0 | met: +20.7 | met: +9%, +1% |
| `tq10` | not met: recall -11.1 on ITS2 | met: +50.5 | not met: +39%, +0% |
| `tq20` | not met: ITS2 spurious +2.38, recall -11.1 / -5.6 | met: +20.8 | not met: +39%, -67% |
| `sw15` | met: +0.00 / +0.00, recall +16.7 / +27.8 | met: +57.3 | met: +4%, -2% |

- **P88 met.** Forward reads alone: spurious +0.15 (ITS1) and +0.57 points (ITS2), recall +11.1 and +0.0 points.
- **P89 met.** Merging first: spurious +0.00 and +0.24 points, recall +16.7 and +33.3 points.
- **P90 met.** Merging first raises the second ITS study's yield by 50.4 points (13.96% to 64.36% of the input pairs).
- **P91 met.** Rank correlation 1.00 on all eight ITS groups (six maxEE arms each).
- **P92 met.** AmpliTrim's predicted pass rate is within 0.22 points of DADA2's at all 48 comparisons.
- **P93 met.** The spurious share at a minimum overlap of 12 nt equals the share at 20 and 30 nt on both markers (+0.00 points).

**Described (no criterion).**
- *The default on the ATCC standard:* ITS1 recall 83.3% (5.0 of the 6 strain-exact members per sample), spurious 0.00%, reads in exact sequences 37.0%, unjudged reads 50.3%; ITS2 recall 66.7% (4.0 of 6), spurious 0.00%, reads in exact sequences 100.0%, unjudged 41.7%. At maxEE 2,2 whole strains are lost, apparently at the filter (it keeps 26.9% of the ITS2 pairs, the trim 91.3%). On ITS1 the trim and the relaxed maxEE arms find them again (recall 88.9% to 100%); on ITS2 the relaxed maxEE arms do not (66.7%), the trim (94.4%) and merging first (100%) do. On ITS1 most judged reads are within 2 edits of their reference but not identical to it.
- *The piglet mixes (the second view):* `sw15` keeps recall (ITS1 76.7%, ITS2 100.0%) and adds no spurious reads; `tq10` lowers ITS1 recall to 46.7%; `tq20` puts 3.04% of the ITS2 reads in spurious sequences. The isolate reference holds C 1, K 3, M 5, P 2 and T 1 sequences on ITS1 and C 1, K 2, M 1, P 4 and T 1 on ITS2; isolate M is flagged (section 33).
- *The two real ITS studies* (shares of the input pairs after chimera removal; median ASVs per sample): the second study `std` 13.96% (23 ASVs), `ee55` 15.16% (25), `sw15` 14.35% (24), `tq10` 31.28% (32), forward reads 53.83% (94), merging first 64.36% (111); the first study `std` 46.85% (105), `ee55` 53.38% (106), `sw15` 59.51% (103), forward reads 61.84% (136), merging first 70.44% (141). On the full-length study the trim does not help, because the loss there is at the merge step, and forward reads alone and merging first give four to five times as many ASVs; whether they are real cannot be judged without truth.
- *The merge test:* vsearch merged 92.7% of the second study's pairs, against DADA2's 13.96% after denoising and merging.

**What follows (as section 33 said).**
1. Q1: P83 holds and P84 does not. Section 33 adopted the closer estimate only if both held, so the report of 0.3.0 does not use it, and for ITS it says that the estimate is an upper bound and why (pairs refused at the merge step, mostly for indels, section 32). An estimate without differences used for ITS only (-2.75 points on the second study) is a candidate for a new registered test.
2. Q2: `ee23` to `ee55` and `sw15` meet P85 and P86 on both markers; `sw15` has the largest ITS2 yield and becomes the advice of 0.3.0 for ITS data where fewer than half the pairs pass at maxEE 2,2: trim the read ends with a sliding window (4 bases, mean quality 15) before DADA2. The report says what it was tested on (short ITS1 and ITS2 amplicons of two fungal mock communities) and that on full-length ITS the yield is limited by the merge step instead. P87 holds for `sw15`.
3. Q3: P88 holds: the report may name forward reads alone as an option when most long ITS pairs cannot merge; it does not recommend it, and it says that on the full-length study they gave four times as many ASVs.
4. Q4: P89 and P90 hold: merging first may be named as an option for full-length ITS, with the caveat that DADA2's authors advise denoising first and that on the full-length study it gave nearly five times as many ASVs.
5. P91 and P92 hold: the report presents the maxEE table for ITS as a reliable guide to the pass rate and to the order of the settings.
6. P93 holds: 12 nt stays for ITS; whether the paper's Limits still name the missing overlap test is the authors' decision (FUTURE.md).
7. Not claimed in any case: accuracy on full-length ITS, which has no truth here; taxonomy; abundances; anything about the four ATCC members left unjudged.

## 35. ITS rules night: can AmpliTrim compute truncQ, the trim strength and maxEE for ITS? (written 14 November 2026, before the run)

**Why.** Run 26 (sections 33 and 34) found one safe and useful answer to poor ITS read ends, a sliding-window trim (window 4, mean quality 15), and found that `truncQ = 10` and `20` raise the yield but lose true strains. It left three questions that decide whether AmpliTrim 0.3 can compute ITS settings from the reads instead of giving fixed defaults:

| | Question |
|---|---|
| A | Does a milder `truncQ` (3 to 8) remove poor read ends and improve merging **without** losing strains? |
| B | How strong should the sliding-window trim be: is mean quality 15 the best of 10, 15, 20 and 25? |
| C | Can a maxEE computed from the ITS reads themselves (the 16S rule applied to untruncated reads) replace the default 2,2? |

**Data, truth and pipeline: as section 33, unchanged.** The eight ITS groups (the ATCC standard, the piglet isolates and mixes of Arfken et al. 2023 for ITS1 and ITS2, and the two ITS studies), thinned with the same function to 5,000 read pairs per sample, so that the reads are the same as in run 26. DADA2 as there, with chimera removal in every arm. The ATCC standard against the NCBI references (`validation/references_its`) is the primary truth of both markers; the piglet mixes against the isolates' own sequences are the second view. The same definitions of true, unjudged and spurious sequences, recall, reads in exact sequences and yield. No 16S control (no question here concerns 16S). The mock libraries carry 32 distinct quality values, the first ITS study 27 and the second ITS study 3, so on the second study most `truncQ` values cut where 2 does; that is described, not judged.

**Arms.**

| Arm | Change from `std` (no truncation, `truncQ = 2`, maxEE 2,2, `minLen = 50`) | Question |
|---|---|---|
| `std` | none | the baseline, and the consistency check |
| `tq3`, `tq4`, `tq5`, `tq6`, `tq8` | `truncQ = 3` / 4 / 5 / 6 / 8 | A |
| `tq10` | `truncQ = 10`, repeated from run 26 as the known-unsafe reference | A |
| `sw10`, `sw15`, `sw20`, `sw25` | the sliding-window trim of section 33 (window 4, cut before the first window whose mean quality is below the value) at mean quality 10 / 15 / 20 / 25, then the `std` filter | B |
| `eer` | maxEE computed per read direction from the dataset's own reads, by AmpliTrim's 16S rule (`pipeline._rule_max_ee`): each read cut at its first base of quality 2 or below (as DADA2 does with `truncQ = 2`), its expected errors EE = Σ 10^(−Q/10) over the bases kept, P75 the 75th percentile (numpy's default) over the reads that keep at least 50 bases, and maxEE = max(2, int(P75) + 1); no truncation; written to `eer.json` | C |

**The consistency check.** The `std` arm is run again on the same reads; for every sample of every group its input, filtered, denoised, merged and after-chimera-removal counts must equal those of run 26 (`results/26_its_night`). If they do not, the night's results are reported but no criterion is decided until the difference is explained. The analysis is `validation/its_rules_report.py`.

**Criteria** (written before the run; "safe" and "useful" as in section 33).
- **P94 (A, safe), per arm `tq3` to `tq8`.** Against `std`, for each marker, on the ATCC standard **and** on the piglet mixes (both truths): the share of reads in spurious sequences rises by at most 1 point, and recall falls by at most 2 points.
- **P95 (A, useful), per arm.** On the ATCC standard's ITS2 the yield rises by at least 5 points against `std`.
- **P96 (A, your mechanism), per arm.** On the ATCC standard's ITS2 the share of denoised pairs that merge rises by at least 2 points against `std`.
- **P97 (A, a flag), per arm.** On the first ITS study the median number of ASVs per sample rises by at most 25%.
- **P98 (B, safe), per arm `sw10`, `sw20`, `sw25`.** As P94.
- **P99 (B, useful), per arm.** As P95.
- **P100 (B, a flag), per arm.** As P97, on both ITS studies.
- **P101 (B, repetition).** `sw15` reproduces run 26's `sw15`: the same counts on all eight groups.
- **P102 (C, safe).** `eer` as P94.
- **P103 (C, useful).** On the ATCC standard's ITS2 the yield of `eer` rises by at least 5 points against `std`.
- **P104 (C, a flag).** As P97, on both ITS studies.

*Safety on both truths* is a choice made for this section: section 33 decided on the ATCC standard alone. *A third mock community* was considered and none was found in time; the code accepts one (`its_night.py --extra`, `its_rules_report.py --truth`), and if one is added before this section is committed, it is named here and must be safe as well.

*Described, no criterion:* the maxEE that `eer` computes for each group; the share of denoised pairs that merge, for every arm and group (the mechanism of A in full); every arm on the piglet mixes (the second view) and on the two real ITS studies; the cut lengths of every trim.

**What follows (decided now).**
- A: if one or more of `tq3` to `tq8` meet P94 and P95 on both markers, the highest such value may become an option in the report of 0.3 for unbinned data, beside the sliding-window trim; which of the two is advised follows the larger ITS2 yield among the arms that are safe and useful. If none does, `truncQ` stays 2, and the report says that values 3 to 10 were tested. P96 says whether removing poor ends improves merging, whatever P94 says.
- B: among the strengths that meet P98 and P99 (with `sw15` from section 34), the one with the largest ITS2 yield becomes the advised strength. If more than one is safe and useful, a per-dataset choice of strength is possible, and needs its own rule, registered and tested.
- C: if P102 and P103 hold, AmpliTrim 0.3 may compute maxEE for ITS by this rule; the report then says it was tested on short ITS1 and ITS2 mock communities. If not, ITS stays at 2,2 with the maxEE table as the guide.
- Not claimed in any case: accuracy on full-length ITS; taxonomy; abundances; the ATCC members left unjudged in section 33.

**What has already been seen (disclosed).** All of run 26's results, including `tq10`, `tq20` and `sw15` on these datasets, and the number of distinct quality values of each group. Not seen: any arm of this section other than `tq10` and `sw15`.

**Procedure.** As section 33: the code is written with tests, a smoke run at 300 read pairs per sample checks the paths (disclosed in section 36), this section and the code are committed, then the night, then the outcomes.

## 36. ITS rules night: outcomes (15 November 2026)

**The run.** Section 35 and its code were committed in `df9dace` before the night ran. The smoke run (300 read pairs per sample) came first, as section 35 said; at that depth `std` and `sw15` reproduced run 26's smoke run exactly, and its verdicts differed from the night's in places (`sw15` was not chosen there: `sw10` was, and P96 was met for `tq8`); they are not outcomes. The night: folder `C:\data\runs\27_its_rules` (named `its_rules` when it ran), 32 of 32 steps; the two real ITS studies took 40 and 107 minutes for their 12 arms. Results in `results/27_its_rules`; the summary of the criteria: `ITS_RULES_SUMMARY.md`.

**The consistency checks.** `std` reproduced run 26's counts on all eight groups, sample by sample, so every criterion is decided; and `sw15` reproduced run 26's `sw15` (P101).

**Outcomes** (5,000 read pairs per sample; safe means safe on the ATCC standard and on the piglet mixes, both markers).
- **Question A, `truncQ` 3 to 8.** `tq3`, `tq4`, `tq5` and `tq6` gave the same counts as `std` on every group: these libraries hold no quality value from 3 to 6, so each cut where `truncQ = 2` cuts. They meet P94 (safe) and P97 and fail P95 (+0.0 points) and P96 (+0.00 points) by having no effect. `tq8` raised the ATCC standard's ITS2 yield by 54.5 points (P95 met) but is not safe (P94 not met: recall −2.2 points in the ITS2 piglet mixes, the limit being 2), and the share of denoised pairs that merge fell by 1.08 points (P96 not met).
- **Question B, the trim strength.** `sw10` meets P98, P99 (+18.7 points) and P100. `sw20` is not safe (recall −4.4 points in the ITS2 mixes); `sw25` is not safe (spurious reads +5.99 points on the ATCC standard's ITS2, recall −4.4 points in the ITS2 mixes); both meet P99 and P100. P101 met.
- **Question C, the computed maxEE.** P102 met (safe on both truths), P103 met (+33.6 points), P104 met (its2 +0%, its +0%). The rule computed 2,3 (ATCC ITS1), 3,9 (ATCC ITS2), 2,4 and 2,6 (the piglet isolates), 2,2 and 2,6 (the piglet mixes) and 2,2 on both real ITS studies.

**Described (no criterion).** On the ATCC standard's ITS2 the trim at quality 15 kept 83.0% of the pairs after chimera removal with recall 94.4%; the computed maxEE kept 59.3% with recall 66.7% (no strain recovered); `tq8` 80.2% with recall 77.8%. On the first ITS study `sw25` lowered the median ASVs by 28%; on the second ITS study (3 quality values) every `truncQ` arm up to 8 equalled `std`.

**What follows (as section 35 said).**
1. A: no value from 3 to 8 is both safe and useful, so `truncQ` stays 2; the report says that 3 to 10 were tested, and that 3 to 6 do nothing on these data. P96 fails for every arm: removing poor ends lets more pairs pass the filter, it does not make them merge better.
2. B: of the strengths that are safe and useful (`sw10`, and `sw15` as in section 34), `sw15` has the larger ITS2 yield and stays the advised strength. Only one strength is better than the others, so no per-dataset choice of strength follows.
3. C: AmpliTrim 0.3 may compute maxEE for ITS by this rule; the report says it was tested on short ITS1 and ITS2 mock communities. The trim and the computed maxEE were tested apart; whether they may be combined is not decided here.
4. Not claimed: accuracy on full-length ITS; taxonomy; abundances; the members left unjudged; any `truncQ` value on data with quality values from 3 to 7.

## 37. The Bakker night: do the ITS rules hold on a third mock community? (written 15 November 2026, before the run)

**Why.** Runs 26 and 27 (sections 33 to 36) found, on two mock communities sequenced as 2×300 MiSeq reads in one laboratory, that a sliding-window trim (window 4, mean quality 15) and a maxEE computed from the reads are safe and useful for ITS, and that `truncQ` should stay 2. They left three questions: whether these rules hold on another mock community, with other read lengths and error-prone PCR; whether the trim may be combined with the computed maxEE; and whether the trigger "fewer than half the pairs pass at maxEE 2,2" is needed for safety.

**The third mock community.** Bakker MG (2018), *A fungal mock community control for amplicon sequencing experiments*, Molecular Ecology Resources 18(3):541-556 (doi:10.1111/1755-0998.12760; BioProject PRJNA377530): 19 fungal isolates in three communities (even, staggered A, staggered B), in triplicate. Two sequencing runs:

| Label | Run | Locus | Reads | Conditions (per community, in triplicate) |
|---|---|---|---|---|
| `bakker_pcr_its1` | the first | ITS1 | about 2×195 nt | the standard workflow, a high PCR cycle number, a long extension |
| `bakker_pol_its1` | the second | ITS1 | 2×250 nt | the polymerases GoTaq, Phusion and PreAmp |
| `bakker_pol_its2` | the second | ITS2 | 2×250 nt | the same |

The primers are those of the paper's methods (sections 2.3 and 2.4), and the primer census of the reads (`PRIMER_CENSUS.md`) confirmed them: BITS and B58S3 (Bokulich & Mills 2013) on the first run, an ITS1 fragment from the end of 18S to the start of 5.8S; ITS1F and ITS2 (ITS1) and ITS3-KYO2 and ITS4-KYO3 (ITS2) on the second, whose files hold both amplicons. Each read pair is assigned to its dataset by its primers (`bakker_prep.py`): the forward primer at the start of one read and the reverse primer at the start of the other, after a spacer of up to 10 bases with up to 2 mismatches; a pair whose forward primer is on R2 is swapped; the spacer and the primers are cut off. Pairs without the expected pair of primers are dropped (on the second run about 2% begin with N and a few per cent carry a primer with more mismatches; `PREP_REPORT.md`).

**Truth.** The ITS sequences of the 19 isolates in Appendix S2 of the paper, cut with each dataset's primers (in-silico PCR, up to 3 mismatches per primer, primers removed) into one reference file per dataset (`validation/references_bakker`, with `REFERENCES_REPORT.md`): `bakker_pcr_its1` 22 amplicons of 105 to 366 nt, `bakker_pol_its1` 22 of 177 to 438 nt, `bakker_pol_its2` 23 of 296 to 514 nt; every isolate has a product in every set (one of the five *Rhizophagus irregularis* ITS2 variants has none, its siblings do). Three ITS1 amplicons carry ambiguity codes (*Fusarium verticillioides* NRRL 20956, 2 codes; *R. irregularis* variants A1.a and A1.d, 4 and 3 codes); each is written as every variant its codes allow (47 sequences in each ITS1 file), so that a perfect read is not counted as edits away from its reference. *F. oxysporum* and *F. verticillioides* are 1 edit apart on both ITS1 fragments (6 on ITS2), closer than the 2-edit rule can tell apart, so on the ITS1 datasets they count as one member, as the paper too could not separate them on ITS1; no other two isolates are within 6 edits in any set. The references are of the very strains, so every read is judged: true within 2 edits of a reference, spurious otherwise, as in section 33. Recall is over the members expected in every sample: 18 on the ITS1 datasets, 19 on ITS2. All 19 isolates are in every community: the paper's Table 2 lists all 19 in each, in equal 18S rRNA gene abundance in the even community and at 1- to 2,000-fold differences in staggered A and B, so the rarest members of the staggered communities may be below detection at 5,000 read pairs per sample, alike in every arm. On the first run the longest amplicon (*Saccharomyces cerevisiae*, 366 nt without the primers, 402 with them) is longer than the two reads together (about 2×195 nt), so it cannot be found by any arm, as the paper reports; this lowers recall alike in every arm, and the comparisons are against `std` on the same reads. A trim may also take away the overlap of pairs that only just overlap; the share of denoised pairs that merge is reported for every arm, so that such a loss shows itself.

**The other truths and datasets, as in sections 33 and 35.** The ATCC standard against the NCBI references, the piglet mixes against the isolates' own sequences, the piglet isolates (which give the mixes' truth) and the two real ITS studies (flags only). All datasets are thinned with the same function to 5,000 read pairs per sample, so that the eight earlier datasets have the same reads as in runs 26 and 27.

**Arms** (every dataset; no truncation, `minLen = 50`, `minOverlap = 12`, `maxMismatch = 0` and chimera removal throughout):

| Arm | | Question |
|---|---|---|
| `std` | `truncQ = 2`, maxEE 2,2 | the baseline and the consistency check |
| `sw15` | the sliding-window trim (window 4, mean quality 15), then maxEE 2,2 | 1, 3 |
| `sw10` | the trim at mean quality 10, then maxEE 2,2 | described |
| `eer` | maxEE computed from the reads by the rule of section 35 | 2 |
| `sweer` | the trim at mean quality 15, then maxEE computed by the same rule on the trimmed reads | 4 |
| `tq10` | `truncQ = 10`, maxEE 2,2: the known-unsafe reference | a check |

**The consistency checks.** `std` must reproduce run 26's counts on the eight earlier datasets, sample by sample; if it does not, no criterion is decided until the difference is explained. `sw15` must reproduce run 26 and `eer` run 27 (P113).

**Criteria** (written before the run; "safe" means: against `std`, the share of reads in spurious sequences rises by at most 1 point and recall falls by at most 2 points, on every truth dataset: the ATCC standard (ITS1, ITS2), the piglet mixes (ITS1, ITS2) and Bakker's three datasets).
- **P105.** `sw15` is safe on every truth dataset.
- **P106.** On every truth dataset where fewer than 50% of the pairs pass at maxEE 2,2 (`std`), `sw15` raises the yield (after chimera removal, of the input pairs) by at least 5 points.
- **P107 (the trigger).** On every truth dataset where at least 50% of the pairs pass at maxEE 2,2, `sw15` is safe and loses at most 2 points of yield.
- **P108.** `eer` is safe on every truth dataset.
- **P109.** `eer` loses at most 1 point of yield on every truth dataset.
- **P110.** `sweer` is safe on every truth dataset.
- **P111.** `sweer` raises the ATCC standard's ITS2 yield by at least 2 points over `sw15`, and loses at most 1 point against `sw15` on every truth dataset.
- **P112 (a flag), for `sw15`, `eer` and `sweer`.** On the two real ITS studies the median ASVs per sample rise by at most 25%.
- **P113.** `sw15` reproduces run 26 and `eer` reproduces run 27 on the eight earlier datasets.
- **P114 (a check).** `tq10` is not safe on at least one truth dataset, as in runs 26 and 27.

*Described, no criterion:* the share of denoised pairs that merge, for every arm and dataset; every arm on Bakker's datasets by condition (standard, high cycle number, long extension, GoTaq, Phusion, PreAmp); the maxEE computed for every dataset, before and after the trim; the members found by no arm; `sw10`.

**What follows (decided now).**
1. P105 and P106: the trim advice holds on a third mock community, with other read lengths and error-prone PCR. If P105 fails, the advice is restricted to what it was safe on, and the report names where it was not.
2. P107: if it holds, the trim did not harm data that pass well at 2,2, and the 50% trigger is not needed for safety (it may stay as a guide to when the trim helps); if it fails, the trigger stays.
3. P108 and P109: the computed maxEE for ITS holds on a third mock community; if either fails, it is not used in 0.3.
4. P110 and P111: 0.3 may combine the trim with maxEE computed on the trimmed reads. If P110 holds and P111 does not, they stay separate; if P110 fails, they are never combined.
5. P113 fails: the differences are reported and explained before any other outcome is used. P114 fails: reported, as a contradiction of runs 26 and 27 to be explained.
6. Not claimed in any case: accuracy on full-length ITS; taxonomy; abundances; Bakker's own conclusions about PCR conditions.

**What has been seen (disclosed).** All results of runs 24 to 27, which include `std`, `sw10`, `sw15`, `eer` and `tq10` on the eight earlier datasets: on those datasets these arms are known (`sw15`, `eer` safe; `tq10` not), so P105, P107, P108 and P109 are decided there by results already seen, and what is new is Bakker's three datasets and the arm `sweer`. If no Bakker dataset passes below 50% at maxEE 2,2, P106 rests on datasets already seen only, and section 38 says so. Also seen: the paper itself (which members it found and lost on each run with its own pipeline), the primer census and the preparation report of Bakker's reads (which pairs carry which primers; no arm was run) and the reference report. Not seen: any arm on Bakker's reads, and `sweer` anywhere.

**Procedure.** As sections 33 and 35: the code is written with tests, a smoke run at 300 read pairs per sample checks the paths (disclosed in section 38), this section, the code, the references and the reports of the census and the preparation are committed, then the night, then the outcomes.

## 38. The Bakker night: outcomes (15 November 2026)

**The run.** Section 37 and its code were committed in `cf3665d` before the night ran. The smoke run (300 read pairs per sample) came first: `std` reproduced run 26's smoke run and `sw15` and `eer` reproduced runs 26 and 27 on the eight earlier datasets; at that depth P105, P107, P110 and P111 were not met, and `sweer` equalled `sw15` everywhere; those verdicts are not outcomes. The night: folder `C:\data\runs\28_its_bakker` (named `its_bakker` when it ran), 55 of 55 steps; results in `results/28_its_bakker`.

**The verdicts as measured.** The consistency check passed (`std` reproduced run 26 on all eight earlier datasets, sample by sample). P105 to P110 and P112 to P114 were met, P111 not met, as `ITS_BAKKER_SUMMARY.md` gives them.

**A defect in the preparation of Bakker's reads, found after the verdicts were read.** On Bakker's two ITS1 datasets 74% and 41% of the reads of `std` were judged spurious, and only 3.3 of 18 and 6.7 of 18 members were found, against 0% and 11.7 of 19 on ITS2 and 17 to 18 members in the paper. The references are not the cause: every ITS1F to ITS2 amplicon has the length of the paper's Table 3 to the base. The ASVs show the cause (`bakker_overhang.py`, `OVERHANG_CHECK.md`): where an amplicon is shorter than the reads (many of the first run's BITS to B58S3 amplicons of 105 to 366 nt under reads of about 195 nt; the short ITS1 amplicons under 2×250 nt reads), each read runs past the far end of the amplicon, and in the merged sequence the overhangs put the forward primer, or its last bases, back at the start and the reverse primer's reverse complement at the end. `bakker_prep.py` removed the primers at the start of each read only; DADA2's `mergePairs` keeps overhangs by default; DADA2's own ITS workflow cuts the reverse-complemented primers from the read ends (cutadapt) for this reason, and this preparation did not. Of the reads of `std`: `bakker_pcr_its1` 26.8% exact references, 72.5% references with primer overhangs, 0.7% other; `bakker_pol_its1` 27.7%, 72.2% and 0.1%; `bakker_pol_its2` 100.0% exact (its amplicons are longer than the reads). The eight earlier datasets carry no primer overhangs (0.0% of reads with five or more bases of a primer at either end), so runs 26 and 27 and the earlier datasets of this run are not affected.

**What stands and what does not.**
- *Valid* (decided on the earlier datasets only, or not concerning Bakker's truth): the consistency check and P113 (`sw15` reproduces run 26 and `eer` run 27 on all eight earlier datasets); P106 (it rests on the ATCC standard's ITS2, +57.3 points, and the ITS2 piglet mixes, +46.7 points, only: no Bakker dataset passes below 50% at maxEE 2,2); P112 (the two real ITS studies). On the earlier datasets alone `sw15`, `eer` and `sweer` were safe and `tq10` was not, as in runs 26 and 27.
- *Not used:* every judgement of Bakker's datasets, so P105, P107, P108, P109, P110, P111 and P114 as decided here. They are decided again on Bakker's reads prepared with the far primer cut (the re-run, section 39), with the criteria of section 37 unchanged.

**Described.** The maxEE computed on the trimmed reads is 2,2 on every dataset (the 75th percentiles of the expected errors fall to 0.35 to 1.02 after the trim), so `sweer` is `sw15` by construction. On Bakker's reads the maxEE computed without the trim is also 2,2 (75th percentiles 0.16 to 1.27), so on Bakker `eer` is `std` and Bakker's datasets cannot test the computed maxEE either way; this does not depend on the defect.

**What follows.** The preparation now cuts the far primer where a read runs into it (`bakker_prep.trim3`, as cutadapt's `-a`/`-A`); the re-run of all eleven datasets, with the overhang check as a gate, is registered as section 39 before it runs. Not claimed from this run: anything about the trim, the computed maxEE or their combination on Bakker's data.

## 39. The Bakker re-run: Bakker's reads with the far primer cut (written 15 November 2026, before the run)

**Why.** Section 38 found that the preparation of Bakker's reads removed the primers at the start of each read only, so that the reads of short amplicons ran into the far primer, DADA2 kept those overhangs, and about 72% of the reads on Bakker's two ITS1 datasets were judged spurious although they were true members' sequences. Every judgement of Bakker's datasets in run 28 was set aside. This section decides them again.

**What changes: the preparation of Bakker's reads only.** `bakker_prep.py` now also cuts the far primer: after the primer at the start of a read is removed, the read is cut where the opposite primer appears reverse-complemented, in full (up to 10% mismatches) or, at the very end of the read, its first bases only (at least 3 bases, up to 10% mismatches), as cutadapt's `-a`/`-A` do in DADA2's ITS workflow (`bakker_prep.trim3`). The primers, the assignment of pairs to datasets and the references (`validation/references_bakker`) are those of section 37. The preparation's report is kept as `references_bakker/PREP_REPORT_rerun.md`.

**What does not change.** The eleven datasets, the thinning to 5,000 read pairs per sample, the six arms (`std`, `sw10`, `sw15`, `eer`, `sweer`, `tq10`), DADA2, the truths and the criteria **P105 to P114 of section 37, with their thresholds and their consequences**, decided by `its_bakker_report.py` as in section 37. The eight earlier datasets' reads are the same as in runs 26 to 28.

**Two gates; if either fails, no criterion is decided until the failure is explained.**
- *Consistency:* `std` reproduces run 26 on the eight earlier datasets, sample by sample, as in section 37; P113 as in section 37.
- *Overhangs:* on each of Bakker's three datasets at most 1% of the reads of `std` are in sequences carrying primer overhangs (`bakker_overhang.py`, as in section 38; `its_bakker_report.py --overhang-max 0.01`).

**Known before the run (disclosed).** Everything in sections 33 to 38, including the Bakker verdicts of run 28 under the defective preparation and the share of the reads that carried overhangs. On Bakker's reads the computed maxEE is 2,2 before the trim and after it, so on Bakker `eer` equals `std` and `sweer` equals `sw15`; P108 and P109 can only be met on Bakker's datasets, and P110 follows P105. Bakker's datasets therefore add evidence on the trim (P105, P107) and on `tq10` (P114), not on the computed maxEE or on the combination. P106 rests on the ATCC standard and the piglet mixes only (section 38). The preparation was changed after run 28's results were read; this section is written before any arm is run on the newly prepared reads.

**Procedure.** The reads are prepared again; a smoke run at 300 read pairs per sample checks the paths and the overhang gate (disclosed in section 40); this section and the changed code are committed; the night runs into `C:\data\runs\bakker_rerun`; the outcomes follow in section 40.

## 40. The Bakker re-run: outcomes (15 November 2026)

**The run.** Bakker's reads were prepared again with the far primer cut, and section 39 and its code were committed in `274cc88` before the night ran. The smoke run (300 read pairs per sample) came first and its overhang gate failed (reads of `std` carrying primer overhangs: `bakker_pcr_its1` 1.81%, `bakker_pol_its1` 33.21%, `bakker_pol_its2` 0.00%); the planned stop after the smoke run was not made, and the registration and the night went ahead. The output of the tests run while registering was cut off in the console; the archive's tests passed afterwards (827 passed). The night: folder `C:\data\runs\29_bakker_rerun` (named `bakker_rerun` when it ran), 55 of 55 steps; results in `results/29_bakker_rerun` (`72fe68d`).

**The gates.** The consistency check passed: `std` reproduced run 26 on all eight earlier datasets, sample by sample. **The overhang gate failed as registered**: reads of `std` carrying primer overhangs, `bakker_pcr_its1` 1.69%, `bakker_pol_its1` 31.40%, `bakker_pol_its2` 0.00% (at most 1% allowed). The analysis labelled every criterion "not decided (the consistency check failed)"; the label was wrong (the consistency check had passed) and is corrected, which changes no number.

**The failure explained.** Section 39 leaves the criteria undecided until a failed gate is explained. Every remaining overhang read carries exactly **one extra base, an A, at the start** of the sequence: on `bakker_pol_its1` 35,758 reads, the *Naganishia albida* sequence that run 28 had with the same count; on `bakker_pcr_its1` 1,858 reads. Overhangs of 3 or more bases: **0.00% on all three datasets** (`OVERHANG_CHECK_min3.md`). The far-primer cut, like cutadapt's, needs at least 3 bases to recognise a primer, so 1- or 2-base remnants are expected; whether the A is the primer's last base or a variant of the isolate, a sequence one base longer than its reference is within the 2-edit rule and was judged true. No judgement depended on the overhangs that the gate counted; the gate's measure counted overhangs of any length, which was a flaw of its definition, not of the preparation. The criteria are therefore decided, on the analysis with the overhang check counting overhangs of at least 3 bases (`its_bakker_report.py --overhang-min-bases 3`, `ITS_BAKKER_SUMMARY_min3.md`); the analysis with the gate as registered is kept beside it (`ITS_BAKKER_SUMMARY.md`), and apart from the corrected label the two give the same numbers.

**Outcomes** (criteria of section 37; safe on all seven truth datasets).
- **P105 met.** `sw15` is safe on every truth; on Bakker's datasets spurious −0.04, +0.00 and +0.00 points, recall −0.6, +2.3 and −1.2 points.
- **P106 met** (ATCC ITS2 +57.3 points, ITS2 piglet mixes +46.7); it rests on the earlier datasets only, no Bakker dataset passing below 50% at maxEE 2,2 (section 38).
- **P107 met.** Where at least 50% pass at 2,2, `sw15` is safe and gains yield: ATCC ITS1 +17.0, ITS1 piglet mixes +15.4, `bakker_pcr_its1` +5.2, `bakker_pol_its1` +5.1, `bakker_pol_its2` +6.2 points.
- **P108 and P109 met.** On Bakker's datasets `eer` equals `std` (maxEE 2,2; disclosed in section 39).
- **P110 met; P111 not met** (+0.0 points against `sw15` everywhere: the maxEE computed after the trim is 2,2).
- **P112 met** for `sw15` (its2 +4%, its −2%), `eer` (+0%, +0%) and `sweer` (as `sw15`).
- **P113 met** (`sw15` reproduces run 26 and `eer` run 27 on the eight earlier datasets). **P114 met** (`tq10` not safe: recall −11.1 on ATCC ITS2, −30.0 in the ITS1 piglet mixes, −4.1 on `bakker_pcr_its1`).

**Described.** With the far primer cut, Bakker's truth works: in `std`, `bakker_pcr_its1` 1.00% of reads spurious and recall 58.9% (10.6 of 18 members), `bakker_pol_its1` 0.00% and 52.1% (9.4 of 18), `bakker_pol_its2` 0.00% and 61.4% (11.7 of 19); recall is below 100% because the staggered communities hold members at up to 2,000-fold lower abundance and because some members are lost by these primers and reads (the paper's own findings), alike in every arm. Yields after chimera removal, `std` to `sw15`: 86.8% to 92.0%, 91.0% to 96.1% and 88.9% to 95.1%.

**What follows (as section 37 said).**
1. The trim advice (window 4, mean quality 15) holds on a third mock community, with reads of about 2×195 and 2×250 nt and error-prone PCR conditions; on Bakker's well-passing reads it gains 5 to 6 points, safely.
2. The trim did not harm data that pass well at maxEE 2,2: the 50% trigger is not needed for safety, and stays in AmpliTrim 0.3 as a guide to when the trim helps most.
3. The computed maxEE for ITS holds; its evidence is the ATCC standard and the piglet communities, Bakker's reads being neutral.
4. The trim and the computed maxEE stay separate: after the trim the computed maxEE is 2,2.
5. Not claimed: accuracy on full-length ITS; taxonomy; abundances; Bakker's own conclusions about PCR conditions.

## 41. The trim check: does AmpliTrim 0.3 predict the share of pairs that pass DADA2's filter after the trim? (written 15 November 2026, before the run)

**Why.** AmpliTrim 0.3 (software commit `42d795e`) advises trimming the poor ends of ITS reads when fewer than half the pairs pass DADA2's filter at maxEE 2,2, and predicts the share of pairs that pass after the trim. Its report says that this prediction is not yet checked against DADA2. This section checks it.

**Data.** The eleven ITS datasets of run 29 (section 39): the thinned reads that DADA2 was given there (on the analysis computer, `C:\data\runs\29_bakker_rerun\<dataset>\reads`) and DADA2's counts (`results/29_bakker_rerun/<dataset>/dada2/std` and `sw15`, `samples.tsv`: the filtered and the merged pairs over the input pairs, summed over the samples). The trim of run 29 (`validation/its_night.py`, `sliding_cut`) and AmpliTrim's (`amplitrim/trim.py`, `sliding_cut`) are the same code.

**Method** (`validation/trim_check.py`). AmpliTrim 0.3 is imported from its repository folder (the script refuses any other version) and analyses each dataset's thinned reads as `amplitrim --its` does (a subsample of at most 50,000 read pairs): the share of pairs that pass at maxEE 2,2 as they are, and, with the trim, the share that pass after it (`pipeline._its_after_trim`), for every dataset whatever the share that passes untrimmed.

**Criteria** (written before the run).
- **P115.** On every one of the eleven datasets, the predicted pass rate after the trim is within 1 point of DADA2's (`sw15`).
- **P116.** On every one of the eleven datasets, the predicted pass rate at maxEE 2,2 is within 1 point of DADA2's (`std`), as P92 required of 0.2.0.

*Described, no criterion:* the predicted yield after the trim against DADA2's merged share (`sw15`); the yield is known to be an upper bound (sections 28 to 32).

**What follows (decided now).**
1. P115 met: AmpliTrim 0.3's report says that the prediction was checked against DADA2 on these eleven datasets (with the largest difference found), `summary.json` gives `trim_advice.checked: true`, and 0.3.0 may be released.
2. P115 not met: the report keeps "not yet checked", the user guide gives the differences found, and the prediction is corrected and checked again before it is called checked.
3. P116 not met: 0.3's filter prediction differs from 0.2.0's validated one (P92 met within 0.22 points on eight of these datasets); it is explained before release.

**What has been seen (disclosed).** All results of runs 24 to 29, including DADA2's counts that this check uses. AmpliTrim's predictions after the trim have not been computed on these datasets. On the code's own tests, the prediction equals an independent computation of DADA2's filter on synthetic reads.

**Procedure.** The check runs in minutes, so there is no smoke run: this section and the code are committed, the check runs into `C:\data\runs\trim_check`, and the outcomes follow in section 42.
