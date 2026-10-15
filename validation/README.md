# Validation toolkit

Manual tools for testing AmpliTrim against real DADA2 runs on public datasets and against simulated ground truth. `pytest` does not
run them. Start with the install in the main README (AmpliTrim and vsearch), then add what these need:

```bat
pip install cutadapt                                    :: removes primers and read-through before AmpliTrim
:: R from https://cran.r-project.org (the tools also look in the standard install folders, `C:\Program Files\R\R-x.y.z`, so a terminal whose PATH lacks R still works), then in R:  install.packages("BiocManager"); BiocManager::install("dada2")
python validation\install_vsearch_windows.py --trust-download     :: Windows only: a checksummed vsearch install
python validation\check_environment.py                            :: shows what is installed and what each route needs
```

The checksum published with vsearch v2.32.0 does not match the file it serves, so the installer stops by default; `--trust-download`
accepts the file (its SHA256 is printed), and `--version 2.31.0` has no checksum to mismatch. Examples below use Windows `cmd`; on
Linux or macOS use `/` in paths and `\` instead of `^` to continue a line.

## Which script is which

Every file of this folder, by what it is for. The scripts import one another by name and start one another by file name, so they stay together in this one folder; this index is the order. Section numbers are those of `TEST_REGISTER.md`. The step-by-step protocol and the details of each tool follow below.

**A. Set up, and check the machine and the build**

| File | What it does |
|---|---|
| `check_environment.py` | Checks Python, AmpliTrim, vsearch, cutadapt, R with DADA2 and FIGARO, and says how to fix what is missing |
| `install_vsearch_windows.py` | Windows: installs vsearch into the active environment, with a checksum |
| `check_build.py` | One place for the build and the run: the build id against `build_manifest.txt`, the machine, what is left to test, fetching and preparing datasets, running everything |
| `example_data.py` | Copies the V4 example reads that ship with DADA2, for a self-test that needs no download |
| `rtools.py` | A helper, not run on its own: finds `Rscript` when R's folder is not on the PATH (ten other scripts import it) |

**B. Get, look at and prepare a dataset**

| File | What it does |
|---|---|
| `select_runs.py` | Shows what a sequencing project contains and picks a small, evenly spaced set of runs |
| `download_ena.py` | Downloads paired FASTQ files from ENA, each checked against ENA's MD5 sum; interrupted transfers resume |
| `describe_reads.py` | Describes paired FASTQ files before analysis: reads, lengths, quality scheme, primers, reversed pairs |
| `longread_describe.py` | The same for long single-end reads: primers and orientation, lengths, expected errors, quality at the ends (section 25) |
| `longread_filter_check.py` | Test L4 (planned; not yet in the register): AmpliTrim's estimate of what a `maxEE` filter keeps against DADA2's `filterAndTrim` on long reads, and how many reads each way of removing primers loses (through `filter_check.R`) |
| `peek_reads.py` | Shows what the reads of some samples look like, for example before and after trimming |
| `real_data.py` | One dataset: removes primers, runs AmpliTrim and, optionally, DADA2 (through `dada2_real.R`) |
| `run_dataset.py` | The whole test protocol on one dataset, collected into one text file to send back (uses `dada2_merge_check.R`) |
| `make_reference.py`, `build_mockA_reference.py` | The expected sequences of a mock community, in the form the accuracy tests use |
| `genome_amplicons.py` | Finds the 16S amplicons of a mock community's strains by searching their genomes with the primers (an in silico PCR) and turns the shares of genomic DNA into the expected share of reads from the copy numbers; writes the form `make_reference.py` writes (long-read truth tests; not yet in the register) |
| `msa1003_genomes.tsv` | The 20 strains of ATCC MSA-1003 with their share of genomic DNA, the NCBI assembly chosen for each, its tier (A exact strain, B another strain) and the reason; read by `genome_amplicons.py` and `night7.py` |
| `night7.py` | Night 7 (planned; not yet in the register): the long-read checks on the PacBio Kinnex-16S data of ATCC MSA-1003 in one run: downloads the genomes, finds the expected amplicons, runs test L4 on three Monomer and three Kinnex files, writes `NIGHT7_SUMMARY.md` (`--check`, `--list`, `--redo STEP` as in `night6.py`) |

**C. Run DADA2, and score what it found**

| File | What it does |
|---|---|
| `compare_methods.py` | Runs DADA2 on the same reads at the settings that different methods recommend (through `dada2_real.R`) |
| `per_sample.py` | Per-sample yield at the recommended settings |
| `max_ee_sweep.py` | DADA2 at one truncation with several `maxEE` values (section 5) |
| `mock_check.py` | Compares the ASVs of a mock community with its expected sequences |
| `classify_asvs.py` | What the ASVs that are not expected sequences are: chimeras, contaminants, error variants |
| `mockA_dilutions.py` | Recovery against input for mock community A |
| `overmerge_real.py` | How many pairs vsearch joins at the wrong place on real reads (section 22) |

**D. The numbered runs, in the order they were run**

The numbers are those of `python validation/organise_data.py --table`, and the same two digits start the folder of each run on the data drive (`05_night3_2026-10-09`). The sections are those of `TEST_REGISTER.md`.

| Run | What | Scripts | Register sections |
|---|---|---|---|
| 01 | Development rounds | `filter_grid_check.py` (test D1), `real_data.py` | 1-8 |
| 02 | Round 4 | `round4.py`, `round4_steps.py` | 9-10 |
| 03 | Night 1 | `night1_v2.py` (not in this archive), `describe_reads.py`, `compare_methods.py` | 11 |
| 04 | Night 2 | `night2.py`, `maxee_arms.py` | 12-13 |
| 05 | Night 3 | `night3.py` | 14-15 |
| 06 | Night 4 (offline) | `maxee_rules.py`, `maxee_frontier.py` | 16-17 |
| 07 | Night 4b | `night4.py` | 18 |
| 08 | ITS check | `its_check.py`, `its_filter_check.R` | 19 |
| 09 | ITS merge check | `its_merge_check.py`, `its_merge_check.R` | 20 |
| 10 | Release check | `release_check.py` | 21 |
| 11 | Soil V4 | `overmerge_real.py`, `overnight.py` | 22 |
| 12 | Merge estimates against DADA2 | `estimate_check.py` | 23 |
| 13 | Night 5 | `nocap_check.py` | 24 |
| 14 | Night 6 | `night6.py`, `longread_describe.py` | 25-26 |
| | The long unattended battery that runs many of the tests above in one command | `overnight.py`, `run_overnight.bat` (the batch file starts it on Windows) | |

`night1_v2.py`, which ran night 1 and prepared its datasets, is **not in this archive**; run 03 can be read from the register but not rerun from here.

**E. Collect the results, and keep the data drive in order**

| File | What it does |
|---|---|
| `collect_results.py` | Copies the reports, protocols, logs and tables of every run folder into the archive's `results/`, and no reads |
| `organise_data.py` | Numbers the run folders and the dataset folders of the data drive in the order of the register, rewrites the folder names in `amplitrim_paths.json`, and finds the files that are stored twice and sets the extra copies aside (never deletes). It plans first and changes nothing without `--apply`; `--check-paths` lists every path of that file as ok or MISSING, before and after the numbering |

**F. Data and the build record**

| File | What it is |
|---|---|
| `PRJNA933120_pairs.csv` | The 26 pairs of runs (volunteer, body site, Illumina run, PacBio run) used by `night6.py` |
| `build_manifest.txt` | The hash of every file the build id covers (`amplitrim/*.py`, `validation/*.py`, `validation/*.R`), read by `check_build.py` |
| `simulated/` | Simulated ground truth: reads with chimeras, indels and repeats, and what scores them (see the section below) |

The R scripts are started by the Python scripts named in the tables (`dada2_real.R` by `real_data.py` and `compare_methods.py`; `dada2_merge_check.R` by `run_dataset.py`; `its_filter_check.R` by `its_check.py`; `its_merge_check.R` by `its_merge_check.py`; `filter_check.R` by `longread_filter_check.py`) and are not run directly.

## Notes on the tools

What the index does not say. The first sentence of each tool is in the index above.

- **`download_ena.py`**: takes a Run Selector CSV, a list, or accessions; header lines such as `acc` are ignored. An ENA server error (HTTP 5xx) is retried for about 90 seconds before it is reported as ENA's, not your connection's. With `--dry-run` a run that cannot be listed does not hide the others. A saved ENA error page is recognised as a web page and named as such. `--max-runs N` takes N runs spread evenly over the list (a pilot, or a small set), and ENA's own run table (its `filereport`, with `run_accession` as the first column) works as the list and, for `--runinfo`, as the sample table. Needs no SRA Toolkit. It takes **paired** runs only: `night6.py` fetches single-end runs with its own small function.
- **`describe_reads.py`**: also gives the **inferred start** of R1 and R2 (a primer shows up as a conserved string with ambiguity letters) and, per sample, which primer the reads start with and after how many bases (barcodes). `--primer SEQ` searches for a primer, `--lengths` gives the length distribution, `--positions N` the base composition.
- **`real_data.py`** (with `dada2_real.R`): removes primers and read-through with cutadapt (optionally sorting mixed-orientation pairs; never quality-trimming), runs AmpliTrim and, optionally, real DADA2 (recommended settings and `--compare F,R`), and appends a row to `results.tsv`. Writes DADA2's per-sample counts and its ASVs. Stops early on a placeholder primer; records a DADA2 failure and carries on.
- **`per_sample.py`**: pooled estimates can hide weak samples. With `--runinfo SraRunInfo.csv` it adds sample names and groups; with `--dada2-track` it adds DADA2's per-sample result and rank correlations.
- **`max_ee_sweep.py`**: with a mock community's expected sequences it scores each result: retention, expected sequences found, share of reads in identical ASVs, extra ASVs.
- **`mock_check.py`**: classifies each ASV as identical, near or none, and says which expected sequences were recovered.
- **`run_dataset.py`**: all of the real-data steps for one dataset in one command, collected into one text file; it runs `dada2_merge_check.R`, which merges the same denoised pairs at 0-8 allowed mismatches and classifies why pairs fail to merge.

**Simulated ground truth** (`simulated/`)

| Script | Purpose |
|---|---|
| `simulate_hard.py` | Simulates reads with chimeras, indels and a perfect repeat across the overlap; every read is labelled with its true template and kind. |
| `eval_hard.py` | Scores every candidate against that truth, per read kind. Imports no AmpliTrim code. |
| `truth_to_tsv.py`, `dada2_check.R` | Convert the truth to a table, and run DADA2 at given truncation lengths and score the ASVs against it. |
| `overlap_sweep.py` | Compares vsearch (strict and tolerant) with DADA2 merge rates across ten overlap regimes. |

## 1. Self-test without downloading anything

```bat
python validation\example_data.py C:\data\example
python validation\real_data.py --local C:\data\example --name dada2_example --out C:\data\results --dada2 --compare 240,160
```

Expected row in `C:\data\results\results.tsv` (identical on Linux with DADA2 1.30.0 and on Windows 11 with DADA2 1.40.0):

| column | value |
|---|---|
| pairs / read_f / read_r | 3000 / 250 / 250 |
| amplicon_median, merged_at_full_length | 253, 0.643 |
| truncLenF / truncLenR, estimated_yield | 151 / 118, 0.8743 |
| dada2_merged, dada2_nonchim, dada2_asvs | 0.88567, 0.88167, 12 |
| compare_dada2_merged (240,160), compare_asvs | 0.62533, 8 |

## 2. Testing a dataset

**Choose.** Illumina paired-end amplicon reads, preferably MiSeq, at most about 40 runs and 1.5 GB, from a paper that states the primers
and region. In the NCBI SRA search (`16S V3-V4 MiSeq 2x300`; filter Illumina, paired, amplicon) open the BioProject, then
Send to > Run Selector > Accession List. Studies deposited directly at ENA (accessions starting `PRJEB`) are held as the submitters'
original files, so they are more likely to have intact quality scores (a heuristic, not a rule). Also usable: the MiSeq SOP mock community,
and mock communities from Zymo or mockrobiota if you have the expected 16S sequences.

**Pilot first.** Many public runs have flattened quality scores (two distinct values), on which DADA2 cannot learn an error model. Download two
runs and describe them before fetching the rest:

```bat
python validation\download_ena.py pilot.txt --out C:\data\NAME\pilot
python validation\describe_reads.py C:\data\NAME\pilot --lengths
```

If `Qvals` is 2 or fewer, DADA2 cannot be tested on this study (AmpliTrim alone still can); choose another. Then download everything
(`--dry-run` lists the sizes first; run it again to resume or finish):

```bat
python validation\download_ena.py runs.txt --out C:\data\NAME\raw --dry-run
python validation\download_ena.py runs.txt --out C:\data\NAME\raw --retries 20
```

With `fastq-dump` instead, use `--split-3` (not `--split-files`), so that R1 and R2 stay in step.

**Run everything with one command** (replace the primers with the ones you found; `FWD` and `REV` are refused):

```bat
python validation\run_dataset.py NAME C:\data\NAME\raw ^
       --primer-f GTGYCAGCMGCCGCGGTAA --primer-r GGACTACNVGGGTWTCTAAT --error-rate 0.15 --compare 240,160 --out C:\data\results
```

`python validation\run_dataset.py --version` prints the build ID of the code you are running (the outcome file starts with it), which tells you at once whether two
copies of the toolkit are the same.

It writes `C:\data\results\NAME_outcome.txt`, even when a step fails, with: the environment; the raw reads described; the trimming, AmpliTrim and DADA2
results; the reads described again after trimming; AmpliTrim's `summary.json`; the per-sample table (16S with 6 or more samples); and, if you ask,
the merge check (`--merge-check`) and the mock check (`--reference FASTA`). Options: `--marker its` (variable-length amplicons; leave out `--compare`),
`--mixed-orientation`, `--pairs N` (pairs per sample given to DADA2 with 6 or more samples, default 5000; use 100000 for low-abundance strains), `--min-reads N` (leave out controls and failed samples with fewer than N pairs), `--only NAMES` (analyse only the named samples: comma-separated, or a file with one name per line; a name that does not exist stops the run),
`--runinfo CSV`, `--explore-coverage` (show what covering fewer of the longest amplicons would allow), `--allow-unmergeable FRAC`, `--explore-max-ee` (also find the best truncation for several `maxEE` values and show what each keeps, the expected errors of the kept reads, and where retention levels off), `--max-ee F,R` (one number, or different forward and reverse limits as in `maxEE = c(2, 5)`; the paper behind a study often uses two). `--compare` is the setting you would normally use: about `240,160` for 2x250 and `270,210` for 2x300.

**Read the output.**

| What you see | What to do |
|---|---|
| `Qvals` of 2 or fewer | DADA2 cannot learn an error model: the run reports "DADA2 not run"; AmpliTrim's result rests on the same flattened quality |
| `Qvals` 3-5 | binned scores (NovaSeq-like): continue, with care |
| a primer at offset 0 (or a few bases in, with a barcode) | pass it as `--primer-f` / `--primer-r` |
| ambiguity letters (Y, R, N, D...) in the "What the reads start with" line | that string is the primer. **The paper can be wrong**: one study named 338F/806R while the reads carried 563F/802R |
| the reverse complement of the other primer in most reads | read-through: the primers option removes it; add `--error-rate 0.15` (cutadapt's default allows only one error in a 15-18 nt primer) |
| `reversed` above 5% | add `--mixed-orientation` (needs the primers) |
| after trimming, primers or reverse complements above about 0% | the trimming did not work; do not trust the rest |
| DADA2 merged far below the estimate | expected for diverse samples at modest depth; run `--merge-check` to see how many pairs fail only by 1-2 mismatches |

Common primers (use the paper's, and confirm them in the output): 341F `CCTACGGGNGGCWGCAG`, 805R `GACTACHVGGGTATCTAATCC`, 515F `GTGYCAGCMGCCGCGGTAA`,
806R `GGACTACNVGGGTWTCTAAT`, 926R `CCGYCAATTYMTTTRAGTTT`, ITS3 `GCATCGATGAAGAACGCAGC`, ITS4 `TCCTCCGCTTATTGATATGC`, ITS1F `CTTGGTCATTTAGAGGAAGTAA`,
ITS2 `GCTGCGTTCTTCATCGATGC`, ITS1F-KYO1 `CTHGGTCATTTAGAGGAASTAA`, ITS2-KYO2 `TTYRCTRCGTTCTTCATC`. A primer found more than 20 bases into the read is a conserved site of the gene (563F lies 29 bases after 515F), not a primer.

**Mock communities.** The MiSeq SOP "Mock" sample (V4, 2x250, primers not in the reads; the DADA2 tutorial page and the mothur MiSeq SOP page link
`MiSeqSOPData.zip`) comes with `HMP_MOCK.v35.fasta`, the expected sequences:

```bat
python validation\run_dataset.py mock C:\data\mock\raw --compare 240,160 --reference C:\data\mock\HMP_MOCK.v35.fasta --out C:\data\results
```

The staggered mock A (PRJEB67797) has 36 runs: ERR12148681-704 are the 24 mock samples, named `KitA-1-<replicate>-<dilution>`, where the
dilution digit is log10 of the input cells (the one failed run, ERR12148698, is replicate 3 at 10^3); ERR12148705-707 are pipeline negative controls and
708-716 PCR negative controls. Below about 10^4 cells over 96% of the reads are contaminants, so analyse the high-input samples (10^8 and 10^9:
ERR12148687, 688, 695, 696, 703, 704) with `--only`, and group the per-sample table by dilution with `--group-regex`; `--runinfo` accepts the
ENA report TSV (`run_accession`, `sample_alias`, `read_count`):

```bat
python validation\run_dataset.py mockA_high C:\data\mockA\raw --primer-f AGAGTTTGATYMTGGCTCAG --primer-r ATTACCGCGGCTGCTGG --mixed-orientation ^
       --error-rate 0.15 --max-ee 2,6 --compare 279,263 --pairs 100000 --explore-coverage --explore-max-ee ^
       --only ERR12148687,ERR12148688,ERR12148695,ERR12148696,ERR12148703,ERR12148704 ^
       --runinfo C:\data\mockA\runs_info.tsv --group-regex "^\d+-KitA-1-\d-(\d)_" --out C:\data\results3
```

Put only the two `Mock_...` FASTQ files in `raw`. The reference may be longer than the amplicon (end gaps are not counted). If a study's mock has several
dilutions, use the high-input samples: below about 10^6 cells most reads are contaminants. To check ASVs later:
`python validation\mock_check.py C:\data\results\mock\asvs_<F>x<R>.fasta reference.fasta`.

**The steps one by one** (for debugging; `run_dataset.py` runs these):

```bat
python validation\describe_reads.py RAW --lengths --positions 15 --primer FWD_SEQ --primer REV_SEQ
python validation\real_data.py --local RAW --name NAME --out OUT --primer-f FWD_SEQ --primer-r REV_SEQ --error-rate 0.15 --dada2 --dada2-pairs-per-sample 5000 --compare 240,160
python validation\describe_reads.py OUT\NAME\ready --lengths --primer FWD_SEQ --primer REV_SEQ
python validation\per_sample.py OUT\NAME --runinfo SraRunInfo.csv --dada2-track OUT\NAME\dada2_per_sample_<F>x<R>.tsv
Rscript validation\dada2_merge_check.R OUT\NAME\dada2_reads <F> <R> 3
```

The runner also accepts accessions directly (`real_data.py runs.txt --name NAME --out OUT ...`). With `--runinfo`, the CSV needs the columns `Run`,
`SampleName` and `spots` (the SRA RunInfo format); sample names are grouped by cutting the replicate off (`--group-regex`, default before `L<number>`).

## 3. Evidence so far

| Dataset | What was checked | Result |
|---|---|---|
| DADA2 example reads (V4 2x250, 3,000 pairs, primer-free) | filter, merge, estimate | filter within 0.1 point of `filterAndTrim`; estimate 87.4% against 88.6% merged at 151/118 (62.5% at 240/160); identical on Linux and Windows |
| One NovaSeq V4 run (2x250, 3 binned quality values, mixed orientation) | filter, estimate | filter exact; estimate 99.9% against 90.1% merged at 127/139 (88.4% at the original 206/60) |
| Soil V4 study, 36 samples (563F/802R, 240 nt reads covering the whole amplicon; the paper named other primers) | filter, recommendation, estimate, ranking | filter matched DADA2 at two settings (98.6% against 98.2-98.8% per sample; 88.5% against 88.4-91.2% at 200/200); recommendation stable from 9 to 36 samples (123/115 then 120/118, total 238); estimate 98.0% against 64.3% merged by DADA2 at 5,000 pairs per sample (47-90% by sample); a long-overlap setting cost 8.0 points in the estimate and 9.8 in DADA2 (9 samples); per-sample rank correlation +0.67 (strict version +0.58); DADA2 retention fell with archive depth (-0.67), cause not established |
| Two public runs (ITS2 and V3-V4, 2x301) | filter | two-value quality scores, so DADA2 could not run; the ITS filter replica matched `filterAndTrim` exactly (148,989 of 157,472 pairs) |
| A 50-sample V4-V5 study (515F/926R, 2x301, 6-nt barcodes) | trimming, amplicon length | two-value quality scores, so DADA2 could not run; trimming removed barcodes, primers and read-through; median amplicon 373 nt, as 515F-926R predicts |
| Staggered mock community A (15 strains, V1-V3 2x300, 17 of 23 samples sequenced in the reverse orientation; 27F-YM/534R) | trimming, filter, estimate, ranking | trimming removed primers and sorted the orientation; the estimate matched DADA2 closely: 65.8% against 63.9% merged at the recommended 275/242 (DADA2 passed 65.6% at the filter, 65.4% were denoised, so denoising lost almost nothing in a low-diversity community); per sample 51.0-69.5% against 50.9-69.2%, rank correlation +0.77 (n = 23), and the one weak sample was flagged by both. The longest amplicon (505 nt) forces a total of 517 nt, so the quality filter, not denoising, is the main loss (about a third of pairs at maxEE 2) |
| Simulated ground truth (`simulated/`) | chimeras, indels, repeats | recommended 220/227 at about 88.5% usable; the default split proposal matched or beat FIGARO on four simulated sets |

| Mock community A at maxEE 2 and 2,6 (V1-V3 2x300, all 23 samples, including low-input dilutions that are mostly contaminants) | effect of a relaxed reverse limit, with the truncation re-optimised | maxEE 2: AmpliTrim 275/242, estimated 65.8%, DADA2 merged 63.9% (50.8% at the authors' 279/263); maxEE 2,6: 257/260, estimated 80.8%, DADA2 merged 77.5% (68.4% at 279/263, so AmpliTrim's truncation retained 9.1 points more under the same filter); ASV counts 229 and 226, so relaxing did not inflate them; per-sample rank correlation +0.77 and +0.63 (strict estimate +0.58 and +0.21) |
| Mock community A at maxEE 2,6, coverage and maxEE tables (all 23 samples, including low-input dilutions) | the cost of coverage and of a looser maxEE | covering up to 512 nt instead of 505 cost 0.2 points of yield (80.6% against 80.8%); leaving out 2-10% of the longest amplicons gave the same setting at 499 nt with a lower yield (79.8%), because more than 10% of the pairs have amplicons of about 499 nt; at 257/260, maxEE 2,2 passed 60.0% and 2,6 passed 81.9%, with the average kept R2 read carrying 0.83 and 1.49 expected errors |
| Mock community A, 6 high-input samples (10^8 and 10^9 cells; V1-V3 2x300), maxEE 2,6 | clean-sample check, the maxEE exploration | 50 ASVs (226 with all samples); the 505 nt floor unchanged and the median amplicon 484 nt (468 with the low-input samples); estimate 80.7% against DADA2's 78.4% merged (73.9% after chimera removal), and 72.3% for the authors' 279/263 under the same filter; per sample 69.1-82.3% estimated, 69.4-80.6% merged; with each maxEE given its own truncation, `3,5` at 268/249 beat `2,6` at 257/260 (83.7% against 81.4% of pairs, 1.80 against 1.92 expected errors per kept pair); the reads `2,6` adds beyond `2,2` carry 3.31 expected errors in R2 |
| Mock community A, 6 high-input samples, real DADA2 at the setting that beats the authors' (`max_ee_sweep.py`, `3,5` at 268/249 against `2,6` at 257/260) | the "beaten by" finding on real data | pass 83.7% (as predicted) against 81.6%; merged 80.8% against 78.4% (+2.4; AmpliTrim predicted +2.2); after chimera removal 76.1% against 73.9%; top-15 share 82.58% against 82.56%; ASVs at 0.1% or more 34 against 34; tail 0.27% against 0.32%; ASVs 54 against 50, the four extra all below 0.1% of the reads. No reference sequences, so the extra ASVs cannot be classed as errors or variants |
| Simulated community where both reads degrade (120 taxa, 430-510 nt) | the "beaten by" rule against real DADA2 | `2,6` at 260/254 was beaten by `3,5` at 273/241; DADA2 at the three settings matched AmpliTrim's pass rates and merged shares exactly (57.6%/56.9% at the default 2,2, 85.1%/83.9% at 2,6, 87.1%/85.9% at 3,5); 114-115 of 120 taxa found, 100% of reads in ASVs identical to the truth and no extra ASVs at any setting |
| Mock community A, 6 high-input samples, real DADA2 sweep at 257/260 (`max_ee_sweep.py`) | does a looser maxEE change what DADA2 finds? | maxEE 2, 2,3, 2,4, 2,5, 2,6: merged 56.4%, 67.1%, 72.9%, 76.3%, 78.4%; after chimera removal 53.3% to 73.9% (chimeras a constant 5.5-5.7% of merged pairs); ASVs 49, 49, 50, 50, 50 (the paper reports 52 mock ASVs for this mock); top-15 share 83.0% to 82.6%; ASVs at 0.1% or more 35 to 34; tail 0.18% to 0.32%. No reference sequences were available, so strain-level accuracy is untested |
| Test set 1: PRJEB42528 (V1-V2, 27F/338R, MiSeq 2x300, 33 samples, 1.43 M pairs), AmpliTrim against FIGARO and against DADA2 (analysis build e3a0f896e8; comparison with the helper scripts of bd59b86e42) | the first head-to-head on a set not used in development | AmpliTrim at maxEE 2: 198/157, estimated yield 90.6%, DADA2 merged 84.2% (82.0% after chimera removal), 808 ASVs. FIGARO's own recommendation 225/130 at maxEE 2,1: predicted retention 88.1%, DADA2 pass 88.5%, merged 82.2%, 857 ASVs. At FIGARO's maxEE (like for like) AmpliTrim chose 211/144 and DADA2 merged 82.3%: a tie with FIGARO (+0.1 point). Two settings with more overlap (230/200: 77.1%; 260/260: 48.3%) merged less, and AmpliTrim ranked all five settings in the order DADA2 did (rank correlation +1.00) with a near-constant optimistic offset of +6.4 points (range +5.7 to +7.7). The filter pass rate was predicted within 0.2 points (91.3% against DADA2's 91.5% at 196/159) |
| Test set 2: PRJEB58189 spring-wheat ITS (ITS1-Fl2/ITS2, MiSeq 2x250, 37 samples, no truncation, maxEE 2), against the authors' own full-depth DADA2 read counts (Table S1 of Yau et al. 2025) | the filter and merge predictions on ITS, and what thinning does | AmpliTrim predicted 71.4% passing the filter; the authors' DADA2 gave a median of 71.1% (69.9% pooled over 3.07 M reads) and ours, on 5,000 pairs per sample, 71.5%. Survival through denoising was 98.2% at full depth against 79.3% on the thinned reads, so thinning cost about 19 points; merged share of the reads entering DADA2 was 56.4% (authors, full depth) against 42.6% (thinned) and AmpliTrim's yield estimate of 69.3%, which is about 13 points above the full-depth figure (DADA2 merged 82.3% of denoised pairs where AmpliTrim assumes about 97%). Across maxEE 1, 2, 3 and 5 the filter prediction was within 0.1 points of DADA2, but the estimate's gain from a looser maxEE (+30 points) was about three times DADA2's merged gain (+13) on the thinned reads. The forward primer was ITS1-Fl2, not the classic ITS1 that had been assumed: trimming with ITS1 left the last 8 bases of ITS1-Fl2 on R1 |
| A trimming defect found on PRJEB58189 ITS (mixed orientation) and fixed | the primer pass of `real_data.py` | with the right primers (ITS1-Fl2/ITS2) AmpliTrim predicted 74.0% passing the filter and DADA2 gave 74.0% (47.2% of the reads merged on 5,000 pairs per sample, against the authors' 56.4% at full depth). But 0.10% of the pairs (3,126 of 3.09 million, the same in forward and reversed samples; the describer's pooled check had suggested 0.5%, which was sampling noise) still had both primers on and the reads the wrong way round: cutadapt counts a read as trimmed when ANY adapter is found, a 3-base chance match of a 3' adapter is enough, and such pairs passed `--discard-untrimmed` untrimmed and the wrong way round (reproduced with a constructed pair). The primers are now found in a strict first pass (at least 8 bases) that decides which pairs are kept, and the read-through is removed in a lenient second pass (down to 3 bases) that discards nothing; a constructed junk pair is discarded, and a primer missing its first base or a 3-base read-through residue are still handled |
| Test set 2, clean reads: PRJEB58189 ITS after the trimming fix (build ba4625e772), 5,000 pairs per sample, maxEE 1, 2, 3 and 5 and a 200/200 truncation | the pre-set ITS criteria on clean reads | filter pass predicted within 0.1 points at all five settings (74.3% against 74.2% at maxEE 2); rank agreement of the estimated yield with DADA2's merged share +1.00; a 200/200 truncation merged 23.3% against 46.9% at maxEE 2 (-23.6 points). DADA2 merged 39.5%, 46.9%, 50.5% and 53.4% against estimated yields of 54.6%, 72.2%, 80.6% and 83.8%: the offset grew from +13.8 to +30.4 points (mean +22.9) with the pass rate, so a constant offset does not describe thinned ITS reads, while the ratio of merged to estimated yield was 0.63-0.72 (mean 0.65; the authors' full-depth ratio of merged to filtered was 0.81). The ratio is a post-hoc observation to be confirmed at full depth. No chimeras were found at any setting (merged equal to non-chimeric), where the authors' pipeline removed 13.9% |
| Test set 2, full depth: PRJEB58189 ITS, all reads of 8 samples (4 forward, 4 reversed; about 680,000 pairs), maxEE 1, 2 and 5, against the 37-sample 5,000-pair run | does depth explain the yield shortfall, and is it constant across maxEE | filter pass predicted within 0.2 points at all three settings (57.9, 75.4 and 95.8% against 57.7, 75.5 and 95.7%). DADA2 merged 51.7, 63.3 and 73.0% against estimated yields of 56.4, 73.1 and 84.2%: shortfall +4.7, +9.8 and +11.2 points (mean +8.5, against +22.9 on 5,000 pairs per sample), so a constant offset again does not describe it, while the ratio of merged to estimated yield was 0.92, 0.87 and 0.87 (0.63-0.72 on 5,000 pairs). Going from maxEE 1 to 5 the estimate gained 27.8 points and DADA2 merged 21.3 more (77%, against 48% on thinned reads); per sample, maxEE 2 beat maxEE 1 by +11.4 points (+10.7 to +12.2) and maxEE 5 beat it by +20.9 (+19.5 to +22.3), in 8 of 8 samples. ASVs 5315, 5456 and 5387, the share of reads in the 15 commonest ASVs 27.1, 25.4 and 24.3%, ASVs at 0.1% or more 146, 149 and 153: no sign that a looser maxEE inflates the ASV set. Merged over filtered was 0.84 at maxEE 2 (the authors' full-depth figure: 0.81; 0.63 on thinned reads). On the 37 samples, per-sample differences against maxEE 2 were +3.7 (+3.3 to +4.1) for maxEE 3, +6.5 (+5.7 to +7.4) for maxEE 5, -7.4 for maxEE 1 and -23.5 (-25.2 to -21.8) for a 200/200 truncation, each in all 37 samples. Chimera removal took out about 0.1% of the merged reads at every setting, where the authors' pipeline removed 13.9%; the R step does detect recurring chimeras in variable-length reads (6.5% of 7.8% planted in a simulation), so the difference is unexplained |
| ITS maxEE table checked against DADA2's filterAndTrim (one simulated ITS-like sample, no truncLen, truncQ 2, minLen 60) | the new ITS table | at maxEE 2,2, 2,5 and 1,1 the number of pairs kept (6862, 8611, 4721), the mean cut length of the kept reads, their mean expected errors for both reads, and the mean expected errors of the pairs a looser limit adds beyond maxEE 2 (0.2901 forward and 3.0129 reverse at 2,5) were identical to the same quantities computed from DADA2's own filtered FASTQ files (largest difference 0); a plain-Python recount in the tests agrees at four settings. A stricter limit shortens the mean kept length (195 to 218 nt from maxEE 1 to 5 in the simulation) because long reads collect the most errors |
| Mock community A accuracy sweep (6 samples of 10^8 and 10^9 cells, V1-V3 2x300, 257/260, build 7dc831b1e5, reference rebuilt by `build_mockA_reference.py`: 53 sequences, 14 strains) | accuracy against a reference, by maxEE | maxEE 2, 2,3, 2,5, 2,6 and 3,5 merged 56.4, 67.1, 76.3, 78.4 and 79.7% of the pairs while the share of reads in ASVs identical to a reference sequence was 99.76, 99.78, 99.78, 99.79 and 99.79%, near 0.02-0.03%, none 0.18-0.21%; the extra ASVs stayed at 14 and the ASV count at 49-50; the chimera share of merged pairs was 5.5% throughout; 35-36 of 53 variants and 11 of 14 strains were recovered exactly at every setting, so the missing three strains are not a maxEE effect. The reference comes from the authors' own DADA2 run (truncLen 299/280, maxEE 2,6) checked against Sanger sequences, so it is not independent truth |
| Mock A: which variants a 257/260 truncation misses (the maxEE 2,6 ASVs of the six 10^8 and 10^9 cell samples, scored against the full 53-variant reference, build 11061686a2) | what limits recovery | 36 of 53 variants and 11 of 14 strains were recovered exactly, but only 39 of the 53 variants have any reads in these six samples (the authors' counts), and 35 of those 39 were recovered. Of the four present variants that were missed, three (Lactobacillus reuteri, 523 nt, 8 and 11 reads) and one (Enterococcus faecalis, 508 nt, 215 reads, 0.22% of the mock reads) are longer than the 505 nt that 257 + 260 can merge with a 12 nt overlap; the fourth (Staphylococcus epidermidis, 497 nt, 36 reads) is mergeable and just rare. The 0.31% of pairs that the AmpliTrim report had flagged as longer than 505 nt are these strains: the coverage rule (up to 1% of the longest amplicons may be left out) lost two real strains, and --allow-unmergeable 0.005 (a total of 524 nt) would have covered E. faecalis |
| Simulated community with a spread of amplicon lengths (120 taxa, 439-504 nt, noisy reverse tail) | the coverage rule (a taxon can merge only if it is at most total - 12 nt) | real DADA2 at totals 511, 505, 480 and 470: the longest taxon recovered was 499, 493, 468 and 457 nt, matching total - 12 to the nucleotide, and none of the 2 + 7 + 68 + 95 taxa longer than that was found; AmpliTrim's estimates were 72.7% and 73.0% against DADA2's 72.7% and 73.3% merged. Giving up long amplicons did not raise the yield here, because a shorter total also means a smaller overlap |
| Simulated mock community (15 strains at 18%/1.8%/0.18%, 2x300 reads, errors exactly as the qualities say, 65% of pairs passing at maxEE 2) | effect of maxEE on DADA2, and the exploration's estimates | at one truncation (275/242), maxEE from (1,1) to no filter: retention 33.5% to 100%, but 15 of 15 strains found, 100% of reads in ASVs identical to the truth and no extra ASVs in every row; AmpliTrim's exact pass rates matched DADA2's (64.8%, 79.7%, 91.8%, 97.2% at 2,2 / 2,3 / 2,5 / 2,8). At each row's own maxEE, AmpliTrim's pass rates equalled DADA2's exactly (47.9%, 77.5%, 93.7%, 96.2%, 96.3%) and its yields were 0.1-1.8 points below DADA2's merged share; the expected errors of the kept reads (R2: 0.885, 1.283, 1.283 at 2,2 / 2,5 / 5,5) equalled the values computed from the reads DADA2 kept. Idealised: real errors are less tidy, so repeat it on a real mock |

What these checks changed in the tools: primers are inferred from the reads, not taken from the paper; cutadapt needs `--times 2` for reads with a primer
and read-through; `--error-rate 0.15` removes read-through with two errors; and a DADA2 failure no longer discards AmpliTrim's own result.

**Not yet validated:** a real mock community, real ITS reads with usable quality scores, V3-V4 on 2x300 with real quality scores, a tight overlap
(V3-V4 on 2x250, where AmpliTrim may report that no setting can cover the amplicon), and DADA2 at the full depth of a diverse study.

## 4. Simulated ground truth

Needs only vsearch (the DADA2 steps also need R). Use a new folder for the simulated reads: the simulator refuses a folder that already contains files.

```bash
D=/tmp/validation
python validation/simulated/simulate_hard.py $D/reads
amplitrim --input-dir $D/reads --work-dir $D/work --out-dir $D/out          # expect truncLenF 220, truncLenR 227, about 88.5% usable
python validation/simulated/eval_hard.py $D                                  # yield_est equals tool_est for every candidate
python validation/simulated/truth_to_tsv.py $D/reads/ground_truth.json $D/truth.tsv
Rscript validation/simulated/dada2_check.R $D/work/subsample/subsample_R1.fastq.gz \
        $D/work/subsample/subsample_R2.fastq.gz $D/truth.tsv $D/dada2.csv 221,226 225,226

python validation/simulated/overlap_sweep.py engines /tmp/sweep
python validation/simulated/overlap_sweep.py dada2   /tmp/sweep              # needs R + dada2
python validation/simulated/overlap_sweep.py report  /tmp/sweep
```

The simulated rates (indels 1e-4 per base, a 32-nt perfect repeat across the overlap, 5% chimeras) are deliberately exaggerated stress levels.

## 5. Does maxEE affect accuracy? (the sweep)

On a mock community, with the trimmed reads of a `run_dataset.py` run, one truncation and several `maxEE` values (DADA2 takes minutes for a few thousand pairs per sample):

```bat
python validation\max_ee_sweep.py C:\data\results\mockA\ready --trunc 275,242 --pairs 6000 --reference expected.fasta ^
       --max-ee 1,1 2,2 2,3 2,5 2,8 --out C:\data\sweep
```

Read it row by row: if retention rises while the expected sequences found, the share of reads in identical ASVs and the extra-ASV count stay the same,
relaxing `maxEE` costs nothing measurable on that dataset. If extra ASVs or non-identical reads grow with `maxEE`, the filter is protecting
accuracy and the value where they start to grow is the limit. Without `--reference` you still see retention, the ASV count and a reference-free abundance check (the share of reads in the 15 most abundant ASVs, the number of ASVs with at least 0.1% of the reads, and the share in ASVs below 0.1%): a mock of N strains should show about N dominant ASVs and a small tail.

## Scoring against a mock community

`mock_check.py` (also used by `max_ee_sweep.py` and by `compare_methods.py --reference`) matches every ASV to its best reference sequence with vsearch:
**identical** (reproduces an expected sequence), **near** (within a few bases, `--min-id`/`--near-id` 0.99, about 4-5 bases over 470 nt, the rule of the paper
that defined the mock A reference), or **none** (a chimera, a contaminant, or a read that is mostly errors). One strain often has several 16S sequences (operon
copies that differ by 1-10 bases), so name the variants of a strain `>ID|STRAIN|label`: the tools then report variants recovered and also the **strains** with at
least one variant recovered (for example 49 of 53 variants, 14 of 14 strains). A near ASV is either an error or a real variant the reference does not list, and
cannot be told apart without the genomes, so report identical, near and none separately.

`build_mockA_reference.py` rebuilds the expected sequences of the staggered mock A (15 strains, V1-V3, ENA PRJEB67797) from the paper's repository
(Huelpuesch, Rauer et al. 2023, BMC Biol 21:269; https://github.com/LuiseRauer/Decontamination_benchmarking, GPL-3.0, so it is downloaded and not redistributed here):
the authors' mock ASVs are the ASVs that contain the middle of a strain's Sanger read (ambiguous calls match either base) plus those within 4 edits. It gives 53
sequences for 14 strains (the paper reports 52; *S. pneumoniae* was not detected). Checks: it contains the seven ASVs the paper names as added by the edit-distance
step; 99.9% of the reads at 10^9 cells and 99.5% at 10^8 fall in these ASVs, 69% at 10^6 and 3% at 10^4; the scoring tool recovers 53 of 53 on the authors' ASV table.
It is the authors' DADA2 ASVs checked against Sanger sequences, not independent truth, and it favours their own settings (truncLen 299/280, trimLeft 20/17, maxEE 2,6). The full reference holds every variant of any of the 23 mock samples, a third of which are not in a given subset of samples and so cannot be recovered
from it; `--samples REGEX` (with `--min-reads`) limits it to the variants present in the samples you analyse, for example `--samples "^(103|104|111|112|119|120)-KitA"`
for the 10^8 and 10^9 cell samples (39 variants, 13 strains).

`mockA_dilutions.py` uses the whole mock A series (23 samples, 10^9 to 10^2 cells, three replicates each): it runs DADA2 on each dilution separately and scores it
against the variants present at that dilution in the authors' own table, giving for every input level the variants and strains recovered, the share of reads in
identical ASVs and the ASVs outside the reference (mostly contaminants at low input). It needs the ENA table (for the sample names), the authors' ASV table
(`Mock_staggered-A_ASVs.csv`, downloaded from their repository) and the reference from `build_mockA_reference.py`.

**`TEST_REGISTER.md` (in the repository root) is the one list of every test: what it asks, the command, the condition (data and pre-set margins) and the outcome, with the test-stage tables by marker and region
and, separately, the development-stage checks (the DADA2 example reads and the simulated communities in `validation\\simulated\\`) that are not results of the manuscript.**

**Can a maxEE rule match FIGARO's? `validation\maxee_frontier.py`.** For a dataset that `compare_methods.py` was run on, it puts FIGARO's own truncation and maxEE on two axes with AmpliTrim's exact quality table (the share of pairs that pass DADA2's filter, and the mean expected errors per kept read), and shows what several candidate maxEE rules would pick at AmpliTrim's truncation. It needs no DADA2 and no vsearch and takes about a minute per dataset: `python validation\maxee_frontier.py --reads <the ready folder> --compare <the compare folder>`. The decision rule applied to its output is in `TEST_REGISTER.md` (decision 2).

**Round 4: the open questions in one place: `python validation\\check_build.py --round4`.** It runs seven steps (registered in `overnight.py` as tier 4, written in `round4_steps.py`) in three phases:
A, quick analyses of results that are already on disk (the forward-heavy split, the variant that only one tool recovers, the pairs a tool needing reads of one length can use, the lengths of the ASVs of the mocks,
and whether FIGARO's maxEE rule can be reproduced); B, simulations with known truth (`simulated/community_sim.py`: how the yield estimate behaves as the number of taxa grows, and how many merges are at the wrong length
with repeats, chimeras and indels); C, the standard battery on any dataset registered since the last run. There is a break of 90 s between the phases (Enter continues, S skips the next phase, Q stops, P holds), and
**the run continues by itself if nothing is pressed**. Start with `--quick --phases A,B` (about 25 minutes). The criteria were fixed in `TEST_REGISTER.md` (section 9) before any run; `ROUND4_REPORT.md` says what each outcome means for the next change.
`python validation/simulated/community_sim.py OUT --taxa 250 --seed 1` makes one simulated community on its own.

## One place: check_build.py

**To run everything that can run without you, in one go (overnight): `python validation\\check_build.py --overnight`.** It also fetches and prepares V4 (PRJNA574030) and both regions of the inter-laboratory mock (PRJEB34919), repeats a step only when a dataset was added since it last ran, keeps Windows awake, and writes `MORNING_REPORT.md` to read in the morning (`--finish` is the same command). It finds your `ready` folders (by their number of samples: 6, 23, 33 and 37 are mock A high input, mock A, V1-V2 and ITS), the input files and `figaro.py`
(skipping the conda and system folders), saves them in `amplitrim_paths.json`, checks the machine, runs every ready step (finished steps are kept; a coverage result saved before the numbers were kept in a readable form is rerun), and prints
the board. `--todo` prints that board alone: each step done, ready or blocked and why; the decision about the default of `--allow-unmergeable`; the new datasets; the manuscript items and the release items that only you can complete; and the
next commands in order. `--decide` applies the rule for the default (change 1% to 0.2% if, on at least two datasets with expected sequences, it recovers at least one more strain on one and fewer variants or strains on none, and across
all datasets it costs a median of at most 2 points of merged pairs and at most 5 on any) to the saved coverage results, and says what is missing when it cannot decide.


`check_build.py` is the single entry point for everything around the validation. `--init` writes `amplitrim_paths.json` (your folders, edited once). Without options it checks the
build id (with `--manifest FILE`, which files differ from another build: a different id means an extra, edited or missing `.py` or `.R` file in `amplitrim/` or `validation/`), the tools
(Python packages, vsearch, cutadapt, R with DADA2), and the data, and prints the `conda install -n <environment> -c conda-forge ...` line for anything missing (for example rapidfuzz).
`--plan` lists every test and download that is left with what each needs. `--download mockB|v4v5|v3v4` fetches the runs from ENA (paired runs only, `--max-runs` spread evenly) and
describes the reads, `--prepare NAME` removes the primers and registers the reads, `--set-reference NAME=FASTA` registers the expected sequences of a mock community, and `--run` runs
`overnight.py` with the paths of the file (`--quick` first).

`classify_asvs.py` says what the ASVs outside the expected sequences are: contaminants (identical or within 2 edits of an ASV that the authors' table finds in the negative controls, with its taxonomy),
chimeras (vsearch uchime, parents among the ASVs or the expected sequences and different from each other), error variants (within 3 edits of an ASV at least 5 times as abundant, or of an
expected sequence), unlisted variants (97 to 99.9% identical to an expected sequence) or unexplained. These are evidence, not proof.

## Overnight run: every test the data on hand supports

`overnight.py` (started by `run_overnight.bat` on Windows) runs eleven steps and writes one `OVERNIGHT_SUMMARY.md`: a table of steps, a table of verdicts against margins fixed in
the script before the run (`MARGINS`, printed in the summary), then the output of each step. A step whose inputs are missing is skipped with the reason, a step that fails is recorded
and the run goes on, and a step that finished is not repeated if you start it again (`--fresh` repeats everything). `--list` shows what is ready, `--only` and `--skip` choose steps,
and `--quick` runs small slices of every step (about 30 minutes) to check the wiring before the long run. Tier 1 steps decide the conclusions of the manuscript; tiers 2 and 3
strengthen it and chase unexplained results.

| Step | Tier | Question | Needs |
| --- | --- | --- | --- |
| mock_grid | 1 | Is AmpliTrim's truncation as accurate as the authors' (variants and strains recovered, reads in identical ASVs, ASVs outside the reference)? Are variants longer than total - 12 nt never recovered, and E. faecalis and L. reuteri recovered once the total allows? | six high-input mock A samples, authors' table, reference |
| asv_classify | 1 | What are the ASVs outside the expected sequences at each setting of mock_grid, and which classes make up the excess of AmpliTrim's truncation over the authors'? | output of mock_grid |
| coverage_levels | 1 | What does leaving out 1%, 0.5%, 0.2%, 0.1% and 0% of the longest amplicons cost (merged pairs) and buy (variants and strains, where there is a reference; ASVs beyond each limit) on every dataset? The evidence for the default of --allow-unmergeable | any dataset |
| landscape | 1 | Is AmpliTrim's recommendation within 1 point of the best of a DADA2 grid of truncations around it that still covers the amplicons? | mock A or V1-V2 |
| extra_datasets | 1 | The same recommendation, coverage levels and comparison with FIGARO and no truncation on new data given with `--extra NAME=READY\|REFERENCE` | new data |
| mock_dilutions | 1 | The same over the dilution series at 257/260, 265/270 and 279/263 | all 23 mock A samples, ENA table |
| mock_clip250 | 1 | A tight-overlap case with known truth: mock A clipped to 2x250 (AmpliTrim, FIGARO, fixed settings) | six samples, reference; FIGARO optional |
| figaro_mockA, figaro_v1v2 | 1 | Does AmpliTrim equal FIGARO at FIGARO's own maxEE (95% interval of the paired per-sample difference within +/- 1 point)? | FIGARO; mock A or V1-V2 reads |
| stability | 2 | Does the recommendation depend on the subsample size? | any reads |
| per_sample_rank | 2 | Per-sample rank agreement between the estimate and DADA2 (at least 10 samples) | 23 mock A samples; V1-V2 |
| mock_pooled | 2 | Does pooled DADA2 inference change recovery? (`AMPLITRIM_DADA2_POOL`) | six samples, reference |
| v1v2_fulldepth | 2 | Calibration of the yield estimate at full depth on V1-V2 (8 samples) | V1-V2 reads |
| its_chimeras | 3 | Is the 0.1% chimera removal (13.9% in the authors' pipeline) a matter of settings? (`AMPLITRIM_DADA2_CHIMERA`, `AMPLITRIM_DADA2_MINFOLD`) | ITS reads |
| mock_miss_trace | 3 | Where are the two unexplained mock A misses lost: absent from the reads, filtered, or denoised? | 23 mock A samples |

`dada2_real.R` reads three optional environment variables (`AMPLITRIM_DADA2_POOL`, `AMPLITRIM_DADA2_CHIMERA`, `AMPLITRIM_DADA2_MINFOLD`); unset, DADA2's own defaults apply.
Not testable without new data: an independent mock community, V3-V4 2x250 and V4-V5 reads, other platforms. `mock_clip250` emulates short reads and is not a substitute for them.

## Comparing methods (for a paper)

`compare_methods.py` runs real DADA2 on one folder of reads with the setting each method would give a user, and puts everything in one table:
what each method predicted, what AmpliTrim's engine says about that setting (so every row is judged by one estimator), and what DADA2 did
(pairs passing, merged, after chimera removal, ASVs, the abundance spectrum and, with `--reference`, expected sequences found and extra ASVs).

```bat
python validation\compare_methods.py C:\data\results3\mockA_high\ready --out C:\data\compare\mockA_high ^
       --amplitrim-max-ee 2 2,6 3,5 --figaro-script C:\figaro\figaro\figaro.py --setting authors=279,263,2,6 ^
       --setting full=0,0,2 --reference mockA_reference.fasta
```

For ITS (`--marker its`; no truncation, FIGARO cannot be used) AmpliTrim's rows are DADA2's ITS workflow at each `--amplitrim-max-ee` value, each with
AmpliTrim's own estimate, so the calibration is shown across `maxEE` values; a fixed truncation copied from 16S shows what it costs:

```bat
python validation\compare_methods.py C:\data\results4\cryobank_its\ready --out C:\data\compare\cryobank_its ^
       --marker its --amplitrim-max-ee 1 2 3 5 --setting trunc200=200,200,2 --pairs 5000
```

Two comparisons, because a tool's recommendation is a truncation *and* a `maxEE`:
- **End to end**: each tool as a user would run it. FIGARO's own top entry is used, with the `maxEE` it pairs with its trim positions.
- **Like for like**: whenever FIGARO is included, AmpliTrim is also run at FIGARO's own `maxEE`, so that only the choice of truncation differs.

`--only NAMES` uses only some samples (comma-separated, or a file of names), and with `--pairs 0` DADA2 then gets all their reads, which is the
way to ask whether thinning to a few thousand pairs per sample cost retention (it can in diverse samples, where DADA2 cannot confirm sequences it sees once).
`peek_reads.py FOLDER [FOLDER ...] --sample NAME` shows the commonest first bases, the lengths and a few example reads of a sample, for example before
and after trimming, to see whether a primer was left behind or which read starts with which primer.

FIGARO needs reads of one length, so it is given the pairs whose reads have the most common lengths (the share is stated in the output); install it
with `pip install scipy` and the FIGARO files from its repository, and pass `--figaro-script`. With six or more samples the output also gives, for each
method against the baseline, the mean per-sample difference in the share of pairs DADA2 merged with a bootstrap 95% interval over samples.

On the simulated community with a spread of lengths (120 taxa, both reads degrading), at equal `maxEE` 2,4 AmpliTrim and FIGARO chose the same truncation
(265/249, 75.1% of pairs merged); end to end they differed only because FIGARO's own `maxEE` (2,4) differs from AmpliTrim's default (2: 277/237, 55.9%).
FIGARO's predicted retention (77.2%) equalled DADA2's pass rate (77.2%) but it does not model merging (75.1%); AmpliTrim's estimates were within 0.7-1.4
points of DADA2's merged share, and it predicted 0.0% for a setting copied from a 2x250 tutorial (240/160), where DADA2 merged 0.0%.

Suggested design, fixed before looking at results: use datasets the tool was **not** developed on as the test sets, and report the others as
development sets; use the same DADA2 steps for every method; report failures (for example reads that FIGARO cannot take) rather than dropping them.

| Role | Datasets |
|---|---|
| Development (used while building the tool) | the DADA2 example reads, one NovaSeq V4 run, the 36-sample soil V4 study, mock community A (V1-V3 2x300), the simulated communities |
| Test (not used so far) | the fungal ITS study PRJDB18761, a V4 mock community with reference sequences (for example the MiSeq SOP mock), a V3-V4 2x300 study with real quality scores, a tight-overlap study (V3-V4 on 2x250) |

## Later runs: the maxEE rule, the ITS checks, the release check

These came after nights 1 to 3. Each has its criteria written down in `TEST_REGISTER.md` before it was run (sections 16 to 21). The scripts import the package of this repository (the development build); the release check and the ITS checks take `--package`, the folder of the clean repository, when they must run the released one.

| Script | What it answers | Needs | Register |
|---|---|---|---|
| `maxee_rules.py` | Does a formula reproduce FIGARO's `maxEE`? The mean and the percentile rules, offline | the datasets | 16, 17 |
| `night4.py` | Night 4b: does the percentile rule merge as many pairs as relaxing to (2,5), and is it as accurate? | R, DADA2, FIGARO, the datasets | 18 |
| `its_check.py --package ...` | Are the ITS report's filter rows and its forward-only share what DADA2's `filterAndTrim` gives, on the same reads? | R, DADA2, vsearch, the ITS reads, night 3's tables | 19 |
| `its_merge_check.py` | Where does DADA2 lose the ITS pairs after the filter, exactly, and does `trimOverhang = TRUE` change it? | R, DADA2, the ITS reads | 20 |
| `release_check.py --package ...` | Does the released package recommend the settings that nights 3 and 4b ran through DADA2, on all 12 datasets? | vsearch, the datasets, the result folders of nights 3 and 4b | 21 |
| `overmerge_real.py` | Over-merging on real reads: how many pairs vsearch joins at the wrong place when the reads are cut to a candidate truncation and merged at the minimum overlap, against a full-length merge as reference; and the errors in the overlap of each cut. Soil V4 by default | vsearch, the dataset | 22 |
| `overnight.py --only landscape --extra "v4=READY"` | Is AmpliTrim's recommendation now near the best of a 3 x 3 DADA2 grid around it? (soil V4: part 2 of section 22) | R, DADA2, vsearch, the dataset | 22 |
| `estimate_check.py` | Which merge estimate is closer to DADA2, strict or tolerant, across all datasets (part A, from the tables on disk, seconds)? And which number of differences vsearch should allow in the overlap to merge as DADA2 does (part B, `--cutoff`, vsearch on the reads)? | the run folders; vsearch and the datasets for part B | 23 |
| `collect_results.py` | Copies the reports, protocols, logs and result tables of the runs into the archive's `results/` (text only, at most 2 MB each, no reads), with a checksum manifest | the run folders | |

Each of them writes its own report (`*_REPORT.md`) to `--out`, and its tests are in `tests/`. Before running anything, `python validation/check_build.py` confirms that this copy is the build whose manifest it holds, and what your machine has.

