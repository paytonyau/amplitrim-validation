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

