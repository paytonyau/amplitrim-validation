# AmpliTrim: future development

Status of 12 October 2026. **Nothing on this list is implemented, and none of it is part of 0.1.0.** Each item says why it is open and what it needs. Anything that would change a recommendation needs a test whose criteria are written in `TEST_REGISTER.md` before it is run, as in nights 1 to 5.

## Before and at the release (logistics, not development)

- Build and upload to TestPyPI, check how the PyPI page renders, then publish to PyPI.
- Switch on Zenodo, tag `v0.1.0`, and put the DOI and the release date in `CITATION.cff`.
- Write the Bioconda recipe (it needs the checksum of the PyPI archive); after the release, add `pip install amplitrim` and the conda command to the README quick start, and a release badge.
- Finish the JOSS paper (`paper/`): fill the three author markers, verify every reference, and submit after the release DOI exists. At submission the paper, its workflow and its test move back to the software repository.
- Fill the seven "source study: to add" entries of the results site, and add the real AmpliTrim report of each dataset (run `validation/collect_results.py`, then rebuild the site).
- Upload the social preview image on GitHub (Settings, Social preview). The brand kit is still in the first colours; regenerate it if the new theme stays.
- Re-check against the primary pages what is still marked "not re-checked" in `docs/FINDINGS.md`: the DADA2 tutorial's advice on relaxing to (2,5), the overlap rule, and the mention of FIGARO as an alternative.

## Algorithm and validation

1. **A simpler cut.** Take the shortest total that covers the amplicons (the longest minus 0.2%, plus the minimum overlap), then the exact best split, and drop the search over longer totals, the yield estimate and the four-key tie-break. Evidence so far (`TEST_REGISTER.md`, section 24): the floor-only cut merged as many pairs as the released cut on all four datasets where the two differ, and 2.3 points more on wastewater (9 of 9 samples). But the pre-set test of equality (P43) was not met (7 of 11), so this needs a new test with its own criteria: all 12 datasets through DADA2 at the floor-only cut, and the stability of the pick over three sample sets.
2. **The `maxEE` rule.** It is a sensible default that has not been tested on independent data.
   - Explain why inter-laboratory V4-V6 was run by night 4b at `maxEE` 3 and is recommended by the released package at 4,3: one unit of forward limit is worth 5.5 points there, and the rule's limit is a step function of a percentile.
   - Validate on new datasets, including rich communities, because what noisier reads do to relative abundances is not measured.
   - Consider a limit expressed as expected errors per 100 nt, which a fixed number cannot be for reads of different lengths.
   - Test another percentile only on data that were not used to choose 0.75.
3. **Repeat the comparison with the other tool for the released package.** The README's comparison (same `maxEE`, the other tool's cut against AmpliTrim's) was run with the development build; the search code is the same and picked the same cuts on ten datasets, but the picks at the other tool's `maxEE` were not compared.
4. **A better estimate.** The estimate stops before DADA2 denoises the reads and is optimistic (on wastewater it preferred the cut that DADA2 liked less). Ideas considered: model the loss in denoising, a strict merge as the headline, an intermediate vsearch cutoff.

## Features

5. **A long-read mode** (PacBio HiFi full-length 16S first, then the whole rRNA operon, 16S-ITS-23S). No truncation and no merging; `maxEE`, `minLen` and `maxLen` from the length distribution of your reads. It needs a single-read input, its own rule (DADA2's author says that keeping about 40% of the reads is fine and that more low-quality reads may not help, so a rule that targets 75% retention does not carry over), and datasets with expected sequences. Nanopore is not claimed: DADA2's author says it works fine with HiFi reads.
6. **A long-amplicon tail diagnostic in the report.** On D6300 V4-V5 about 1% of the amplicons are 550 nt against a median of 374, and leaving out 1% instead of 0.2% raised the merged share from 57.2% to 89.8%. A report line saying what leaving out 1% would allow, and the coverage setting as an option, would let users decide.
7. **Smaller ideas, each set aside for 0.1.0:** a `maxEE` adjustment for ITS; a forward-heavier split; a minimum-overlap sweep; a per-cell DADA2 flow; Deblur support (it takes a trim length only); a running-time comparison on the released build (not yet timed).

## Decided against for now

A donation button, and analytics on the website.
