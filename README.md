# AmpliTrim: validation archive

This is the research record behind AmpliTrim's validation: the scripts, the criteria written down before each run, the outcomes, the results of every run, and the write-up. **It is not the software.** To use AmpliTrim, install it from <https://github.com/paytonyau/amplitrim>, whose [`docs/validation.md`](https://github.com/paytonyau/amplitrim/blob/main/docs/validation.md) summarises what held and what did not.

| Path | Contents |
|---|---|
| `docs/FINDINGS.md` | The write-up: methods, the 12 datasets, results, limitations |
| `docs/SUPPLEMENT.md` | The four development rounds that came first |
| `docs/validation.md` | The evidence behind each recommendation of the released package: what held, what did not, and the limits. The software's README links here |
| `paper/` | The draft of the JOSS paper, with its tests (`tests/test_paper.py`) and a workflow that builds the PDF. It moves to the software repository for submission |
| `TEST_REGISTER.md` | Every test: the criteria and predictions written before each run, and the outcomes (sections 1 to 23) |
| `validation/` | The scripts that produced every number; `validation/README.md` has the step-by-step protocol and a table of the later scripts |
| `amplitrim/` | The development build that was tested (not the released package; see below) |
| `tests/` | The package's tests and the tests of the validation scripts |
| `results/` | The reports, protocols, logs and result tables of the runs, with `MANIFEST.tsv` (path, size, sha256); no reads |

## Which software was tested

The study ran a **development build** of AmpliTrim, which this folder holds in `amplitrim/`. Nights 1, 2 and 3 each ran the whole set of datasets. Night 1 was the standard battery (build `6b0c25b6e3`), night 2 the maxEE head-to-head (`b8469465fc`), and night 3 maxEE 2,5 and a careful manual cut (`194b5ebf74`). Night 4b (the percentile rule) and the ITS checks came after. **No file of the package (`amplitrim/`) changed between any of these runs;** only validation scripts were added, which is why the build id changed. This folder holds the last build (`a75392d6f8`), and `python validation/check_build.py` confirms that this copy is exactly that build and checks which tools your machine has.

The released package (0.1.0, in the main repository) is **not** identical to it. Compared with this build:

- **Identical:** `__init__.py` and `config.py`. `trimmer.py`, the quality filter and the pass-rate table, only gained two methods (`cut_errors`, `untruncated_forward`); nothing in it was changed. `merge_wrapper.py` and `subsample.py` differ by a comment each, and `ranking.py` lost only the writer of `candidates.csv`.
- **Removed:** the FIGARO cross-check (`figaro_wrapper.py`, its compatibility shim, its options), every command-line option except the three required ones, `--its` (which replaces `--marker its`) and `--verbose`, and the `candidates.csv` output. The validation scripts in this archive use these.
- **Added:** `maxEE` is recommended by a rule. Per read, it is the limit that covers three quarters of the reads at the recommended cut, never below 2 and never above 5. This is the percentile rule of `TEST_REGISTER.md` sections 17 and 18, chosen by the researcher, and it replaces rule R of night 3. Its DADA2 confirmation (night 4b) did not meet its pre-set adoption criteria, so we keep it as a trade-off that we state openly. `summary.json` has `schema_version`, `maxEE_f`, `maxEE_r` and `maxEE_rule`. For ITS, the report shows where the pairs go (`its` in `summary.json`).
- **Changed:** the search and the report (`pipeline.py`, `cli.py`, `report.py`): the FIGARO branches are gone, and the report text and the ITS table are reworded. The per-`maxEE` and per-coverage searches still exist in the code but are only reachable from Python.

DADA2's result depends only on the settings it is given (`truncLenF`, `truncLenR`, `maxEE`, with `minOverlap` and `minLen`). What ties the released package to the runs of this study is whether it **recommends the same settings** as the ones that were run through DADA2. `validation/release_check.py` runs the released package on the 12 datasets and compares them (`TEST_REGISTER.md`, section 21). Its outcome is recorded in the register and in `results/`.

## Results

`results/` holds the text outputs of each run: the reports (`*_REPORT.md`), the protocols (`*_PROTOCOL.md`), the logs, and the result tables and JSON files, with their folder structure. It was filled by `validation/collect_results.py`. The script leaves out the reads, everything that DADA2 and vsearch wrote to work folders, and the ASV FASTA files. It caps each file at 2 MB and replaces a Windows user folder by `C:\Users\<user>`. `results/MANIFEST.tsv` lists every file with its size and sha256, so the record can be checked.

## The data

**The reads are not in this repository.** They are public; the accessions are in Table 1 of `docs/FINDINGS.md`, and `validation/README.md` has the steps to fetch and prepare each dataset. The README of the software shows how to make a small simulated dataset for trying it.

## Reproducing it

You need Python, vsearch, cutadapt, R with DADA2 and, for the comparison, FIGARO. The steps, and the commands for each dataset, are in `validation/README.md`.

## Licence

MIT (see `LICENSE`), as the software. Please cite the software (see the main repository) and, once archived, this record.

**Third-party software.** The scripts here run other programs as separate processes: vsearch (GPLv3 and BSD-2-Clause), cutadapt, R with DADA2, and FIGARO (GPLv3, <https://github.com/Zymo-Research/figaro>) for the comparison. None of their source code is part of this repository, and their own licences govern their use and distribution.
