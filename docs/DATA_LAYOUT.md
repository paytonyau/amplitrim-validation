# Where the data of the validation are, and what is in this repository

The repository holds the scripts, the tests, the register, the written results and the site. It does **not** hold the reads, the DADA2 and vsearch
work folders, the ASV FASTA files or the long-read data: they are large and can be made again. This page says where each part is on the data drive
(`C:\data` in the runs of 2026) and how to rebuild what is missing. Sizes are those of 4 November 2026.

## The data drive

| folder | what it is | size |
|---|---|---|
| `01_...` to `12_...` | the 12 datasets, one folder each, named `NN_<kind>_<region>_<name>_<accession>` (the table: `python validation/organise_data.py --table`); dataset 13 is the dataset of night 6 (PRJNA933120), 14 the MiSeq SOP mock (`miseqsop_mock`, Night 10) and 15 the second ITS dataset (PRJEB108994, Night 10) | 9.7 GB |
| `runs\` | the numbered runs, one or more folders each, named `NN_<what>_<date>` | 11.1 GB, of which 7.9 GB bulk |
| `longread\` | the long-read data of nights 7 to 9b | 4.5 GB |
| `_archive_2026-10\` | deleted on 6 November 2026 (the outputs of the first rounds, superseded by the numbered runs) | - |
| `_release_acceptance\` | the rehearsal and the acceptance run of the PyPI release of 0.1.0 (not validation runs) | small |
| `_scratch\` | test folders and first looks that are not runs of the register | small |
| `backup\`, `_reports\`, `tools\` | copies of `amplitrim_paths.json`, the inventories and plans, and the helper scripts | small |

`amplitrim_paths.json` (the paths of the datasets, their primers and references) is machine specific and is **not** in the repository.

## The datasets

The reads are public and are not in this repository. Re-fetch them with `validation/download_ena.py` (it checks the ENA checksums) and prepare them with the scripts of group B in `validation/README.md`; `python validation/organise_data.py --check-paths` lists which paths of `amplitrim_paths.json` exist.

| Dataset | Accession | Folder on the data drive | Key in `amplitrim_paths.json` |
|---|---|---|---|
| Mock community A (V1-V3) | PRJEB67797 | `01_mock_V1-V3_mockA_PRJEB67797` | `mockA_high`, `mockA_all` |
| Mock community B (V1-V3) | PRJEB67827 | `02_mock_V1-V3_mockB_PRJEB67827` | `mockB` |
| Inter-laboratory mock (V1-V2) | PRJEB34919 | `03_mock_V1-V2_interlab_PRJEB34919` | `mcm_v12` |
| Inter-laboratory mock (V4-V6) | PRJEB34919 | `04_mock_V4-V6_interlab_PRJEB34919` | `mcm_v46` |
| ZymoBIOMICS D6300 (V3-V4) | PRJNA1261676 | `05_mock_V3-V4_zymoD6300_PRJNA1261676` | `d6300_V3-V4` |
| ZymoBIOMICS D6300 (V4-V5) | SRR16961581 to SRR16961583 | `06_mock_V4-V5_zymoD6300_miseq` | `d6300_V4-V5` |
| Human gut (V1-V2) | PRJEB42528 | `07_real_V1-V2_PRJEB42528` | `v1v2` |
| Spring wheat soil (V4) | PRJNA574030 | `08_real_V4_soil_PRJNA574030` | `v4` |
| Spring wheat rhizosphere (ITS) | PRJEB58189 | `09_real_ITS_cryobank_PRJEB58189` | `its` |
| Wastewater (V3-V4) | PRJEB75155 | `10_real_V3-V4_wastewater_PRJEB75155` | `wastewater_V3-V4` |
| Mouse gut with spike-ins (V3-V4) | PRJNA1243306 | `11_real_V3-V4_spikein_PRJNA1243306` | `spikein_V3-V4` |
| Surface swabs (V3-V4) | PRJNA1261676 | `12_real_V3-V4_swabs_PRJNA1261676` | `swabs_V3-V4` |
| Dataset of night 6 (undecided: the Illumina reads have two quality values) | PRJNA933120 | `13_...` | (night 6 only) |
| MiSeq SOP mock community (final run) | none: the mothur MiSeq SOP data (Kozich et al. 2013) | `14_mock_V4_miseqsop` | `miseqsop_mock` |
| Soil fungi, full-length ITS (final run; nine runs by a fixed rule) | PRJEB108994 | `15_real_ITS_fullITS_PRJEB108994` | (given to `night10.py` with `--its-new`) |
| Single isolates of Arfken et al. 2023, ITS1 (follow-up; five isolates in triplicate) | PRJNA693350 | `16_mock_ITS1_arfken_isolates_PRJNA693350` | (given to `its_followup.py`; label `arfken_its1`) |
| Single isolates of Arfken et al. 2023, ITS2 | PRJNA693350 | `17_mock_ITS2_arfken_isolates_PRJNA693350` | (label `arfken_its2`) |
| Mixes MixA, MixB and MixC of Arfken et al. 2023, ITS1 (follow-up) | PRJNA693350 | `18_mock_ITS1_arfken_mixes_PRJNA693350` | (label `arfken_mix_its1`) |
| Mixes MixA, MixB and MixC of Arfken et al. 2023, ITS2 | PRJNA693350 | `19_mock_ITS2_arfken_mixes_PRJNA693350` | (label `arfken_mix_its2`) |
| The ATCC mycobiome standard MSA-1010 of Arfken et al. 2023, ITS1 (overnight run; ten strains in triplicate) | PRJNA693350 | `20_mock_ITS1_arfken_atcc_PRJNA693350` | (label `atcc_its1`; runs in `arfken_atcc_runs.tsv`) |
| The ATCC mycobiome standard MSA-1010 of Arfken et al. 2023, ITS2 | PRJNA693350 | `21_mock_ITS2_arfken_atcc_PRJNA693350` | (label `atcc_its2`) |
| The fungal mock communities of Bakker 2018 (19 isolates; even, staggered A and B), ITS1, the PCR conditions (about 2×195 nt) | PRJNA377530 | `22_mock_ITS1_bakker_pcr_PRJNA377530` | (label `bakker_pcr_its1`, in `prepared/its1/ready`; runs in `bakker_runs.tsv`) |
| The same communities, the polymerases (2×250 nt; ITS1 and ITS2 in the same files, split by primer) | PRJNA377530 | `23_mock_ITSmix_bakker_polymerase_PRJNA377530` | (labels `bakker_pol_its1` and `bakker_pol_its2`, in `prepared/its1/ready` and `prepared/its2/ready`) |

## The runs

The number of a run is its number in `validation/README.md` and in `python validation/organise_data.py --table`: the order of the register, not of the dates.
`organise_data.py` gives each folder its number and rewrites the names in `amplitrim_paths.json` in the same step; it plans first and has an undo.

| run | folders under `runs\` | register |
|---|---|---|
| 01 | `01_development_rounds_2026-10-07` (text outputs only, recovered from the supplement package of 7 October 2026; no reads, no work folders) | 1-8 |
| 02, 10 | no folder under `runs\` (no raw outputs of these runs are kept in this archive; their results are recorded in register sections 9-10 and 21 and in `docs/FINDINGS.md`) | 9-10, 21 |
| 03 | `03_night1_2026-10-09` | 11 |
| 04 | `04_night2_2026-10-09` | 12-13 |
| 05 | `05_night3_2026-10-09` | 14-15 |
| 06 | `06_night4_offline`, `06_night4b_offline`, `06_night4c_offline` | 16-17 |
| 07 | `07_night4_2026-10-10` (night 4b) | 18 |
| 08, 09 | `08_its_check`, `09_its_merge_check` | 19, 20 |
| 11 | `11_overmerge_v4`, `11_v4_landscape` | 22 |
| 12 | `12_estimate_check`, `12_estimate_check_cutoff` | 23 |
| 13 | `13_night5_*`, `13_nocap_check`, `13_nocap_check_mockB`, `13_nocap_check_nights` | 24 |
| 14 | `14_night6_pilot`, `14_night6_PRJNA933120` | 25-26 |
| 15, 16 | `15_longread_describe_sq2_monomer_all96`, `16_longread_describe_sq2_kinnex` | - (`FUTURE.md`) |
| 17 to 20 | `17_night7_longread`, `18_night8_truth`, `19_night9_revio`, `20_night9b_revio` | - (`FUTURE.md`) |
| 21 | `21_release_check_020` (release check of version 0.2.0) | 28 |
| 22 | `22_timing_0.2.0` (run time and memory of 0.2.0) | - |
| 23 | `23_overlap_test_2026-11-04` (the minimum-overlap test: 36 run folders, the batch files, the logs, `overlap_results.*`) | 27 |
| 24 | `24_night10_closing_<date>` (the closing run: the command line on every dataset, the coverage setting, the new datasets, the hold-out overlap test, the grid; `NIGHT10_SUMMARY.md`) | 29-30 |
| 25 | `25_its_followup`, `25_its_followup_arfken`, `25_its_depth_1000`, `25_its_depth_20000`, `25_its_depth`, `25_its_followup_summary` (the ITS follow-up: the pairs that DADA2 refuses, the truth test on the Arfken mixes, AmpliTrim's estimates against DADA2, the depth series; `ITS_FOLLOWUP_SUMMARY.md`) | 31-32 |
| 26 | `26_its_night` (the ITS overnight run: every arm of every dataset, the stricter estimates; `ITS_NIGHT_SUMMARY.md`) | 33-34 |
| 27 | `27_its_rules` (the ITS rules night: truncQ, the trim strengths and the computed maxEE on every ITS dataset; `ITS_RULES_SUMMARY.md`) | 35-36 |
| 28 | `28_its_bakker` (the Bakker night: the trim, the computed maxEE and their combination on three mock communities; `ITS_BAKKER_SUMMARY.md`) | 37-38 |
| 29 | `29_bakker_rerun` (the Bakker re-run: the eleven datasets again, Bakker's reads with the far primer cut; `ITS_BAKKER_SUMMARY.md`, `OVERHANG_CHECK.md`) | 39-40 |
| 30 | `30_trim_check` (the trim check: AmpliTrim 0.3's predicted pass rates against DADA2's on run 29's eleven datasets; `TRIM_CHECK_SUMMARY.md`) | 41-42 |
| 31 | `31_v04_night` (the v0.4 night: twelve ITS datasets with Pauvert's, twenty-two arms, the 16S regression; `V04_SUMMARY.md`) | 43 |

Runs 15 to 20 were planned in `FUTURE.md` (the long-read test plan), not in numbered sections of the register.

## What is in the repository for each run

`results\` holds the text outputs of the runs: the reports, protocols, logs and result tables (`.md .tsv .json .log .txt`, at most 2 MB each), with their folder
structure, and `results\MANIFEST.tsv` (path, size, sha256) so that the record can be checked. `python validation/collect_results.py --runs C:\data\runs --out results --dry-run`
shows what it would copy; it leaves out reads, work folders and ASV FASTA files, and replaces `C:\Users\<name>` by `C:\Users\<user>`.
The ASV files are not kept; the overlap test keeps a fingerprint of every ASV set instead (`overlap_results_asv_sets.tsv`): equal fingerprints mean exactly equal ASVs.

## How to rebuild what is not here

1. The reads: `validation/download_ena.py` and the accessions of `docs/SUPPLEMENT.md`; `python validation/check_environment.py` checks the machine (R, DADA2, vsearch, cutadapt).
2. The prepared datasets (primer-trimmed pairs): the scripts of group B in `validation/README.md`; the primers and references of each dataset are in `amplitrim_paths.json`.
3. A run: the script of its row in the run table, with the command of its register section. The runs that took a 5,000-pair sample per sample say so in their section.
4. A new run gets the next number: `python validation/organise_data.py --base C:\data --next`.

## Keeping it in order

`python validation/inventory_data.py --base C:\data --repo <this repository> --out <report>` lists every folder with its size, its newest file, its bulk against its
results, and flags the run folders without a number, the scratch folders and the runs not yet collected into `results\`. It changes nothing.
