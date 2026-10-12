---
title: 'AmpliTrim: choosing DADA2 truncation lengths and maxEE from the reads themselves'
tags:
  - DADA2
  - amplicon sequencing
  - 16S rRNA
  - ITS
  - microbiome
  - Python
authors:
  - name: Payton (Tung On) Yau
    orcid: 0000-0002-3283-0370
    corresponding: true
    affiliation: 1
affiliations:
  - name: London School of Science & Technology, London, United Kingdom
    index: 1
date: 10 October 2026
bibliography: paper.bib
---

# Summary

DADA2 [@callahan2016dada2] infers amplicon sequence variants from Illumina reads. For paired-end data it asks the user to choose how many bases to keep from each read (`truncLenF`, `truncLenR`) and how many expected errors a read may carry (`maxEE`). These settings decide how many read pairs survive quality filtering, how many of those can still be merged, and whether the longest amplicons can be merged at all. They are usually chosen by eye from quality plots.

AmpliTrim is a Python command-line tool that recommends these settings from the user's own reads. It measures the amplicon lengths, replicates DADA2's quality filter exactly from the quality scores for every forward and reverse split, and tests real merges with vsearch [@rognes2016vsearch]. It writes a report with settings that can be pasted into DADA2, in R or through QIIME 2's DADA2 plugin [@bolyen2019qiime2], a figure of what each choice keeps, and a machine-readable summary for pipelines. A mode for amplicons of very different lengths (ITS) is included. AmpliTrim was tested against real DADA2 runs on 12 public datasets, with the criteria written down before each run.

# Statement of need

Truncating too little leaves low-quality tails that fail the expected-error filter. Truncating too much leaves too little overlap, and merging collapses. The collapse is abrupt: on a mock community, 74 to 76% of pairs merged at the shortest total length that covers the amplicons, and 25% at 40 nt below it. A strict `maxEE` can also cost most of the data silently: on another mock community, 19.9% of the input pairs merged at DADA2's default `maxEE` of 2, against 55.0% at the value AmpliTrim recommends. The usual guidance is visual, from the quality plots of the DADA2 tutorial [@dada2tutorial], with rules of thumb from the author of DADA2 [@dada2issue236; @dada2issue2136]. That is hard to reproduce, hard to justify in a methods section, and hardest for new primer sets, new read lengths, reads whose lengths differ after primer removal, and ITS.

AmpliTrim is for researchers who run DADA2, directly or through QIIME 2, and want settings chosen by a procedure they can describe and repeat, together with what each choice costs.

# State of the field

The DADA2 tutorial [@dada2tutorial] and its author's guidance recommend reading the quality profiles and keeping enough sequence for the reads to overlap. A rule of thumb is a total of the amplicon length plus 25 nt, and the author suggests trying a few `maxEE` values until enough reads pass [@dada2issue236]. A `maxEE` of 2 is described as a good general-purpose value that can be relaxed somewhat [@dada2issue2136]. FIGARO [@weinstein2019figaro] automates the choice of truncation lengths by searching for the best retention under a read-quality model. It expects reads of one length, so after primer removal it could use between 8% and 98% of the pairs in our datasets, and it does not apply to ITS. Read trimmers such as cutadapt [@martin2011cutadapt] and fastp [@chen2018fastp] remove adapters, primers and low-quality ends, but they do not choose DADA2's settings; AmpliTrim assumes that the primers have been removed.

AmpliTrim does not replace FIGARO's truncation search so much as work from different inputs: it uses all read pairs of any lengths, tests merging directly, recommends `maxEE` as well, and has an ITS mode. At the same `maxEE` its truncation merged as many pairs as FIGARO's on our datasets (equivalent, not better), which is the claim it makes.

# Software design

AmpliTrim takes a deterministic subsample of the reads and merges the overlapping pairs to measure the amplicon lengths. The shortest total length that can merge the longest amplicons (all but the longest 0.2%) is computed exactly: the amplicon plus the minimum overlap of 12 nt, which is DADA2's default for `mergePairs`. Totals are tested upward from there. For each total, several forward/reverse splits with a high filter pass rate are given a DADA2-style filter and a vsearch merge test. Among the settings whose estimated yield is within a small tolerance of the best, the one with the shortest overlap is recommended, because it relies least on the noisy ends of the reads.

Two design decisions follow from what can and cannot be exact. The quality filter can be replicated exactly, so it is: its predicted pass rate was within 0.8 points of DADA2's on all 12 datasets. The merge yield cannot, because it is estimated with vsearch before DADA2 denoises the reads. It is therefore reported as an upper bound for comparing settings, and it overshot the merged share by up to 7 points, and by 21 and 32 points on an ITS and a soil dataset. The merge is deliberately tolerant of differences in the overlap (up to 10, vsearch's default), because DADA2 merges denoised reads. A strict merge of the raw reads, with no difference allowed, was farther from DADA2 (12.5 points on average against 6.5) and ranked the settings worse.

`maxEE` is recommended per read. It is the integer part of the 75th percentile of the reads' expected errors at the recommended cut, plus 1, never below 2. There is no upper limit, because the author gives no number; the report notes a limit above 5. Clean reads therefore stay at DADA2's default, which its author calls a good general-purpose value [@dada2issue2136], and the limit is relaxed only when the reads need it, as his advice that relaxing it somewhat is acceptable allows. The percentile is the author's own choice, made after seeing the validation data; they are a trade-off that we state openly, not a result (see below). In ITS mode there is no fixed truncation length, as in DADA2's ITS workflow, and the report shows where the pairs go (too short, containing an N, over `maxEE`, passing but not merging) and what the forward reads alone would keep.

The interface is deliberately small: three required arguments and two options, `--its` and `--verbose`. The outputs are a Markdown report, a figure and a JSON summary with a schema version. The software is tested with simulated reads of known origin, with continuous integration on Linux, macOS and Windows (Python 3.9 and 3.12), and its README shows how to make a small simulated dataset for trying it.

# Research impact statement

The evidence that AmpliTrim does what it claims is its validation. The validation used real DADA2 runs on 12 public datasets (six mock communities and six real datasets, ITS among them). Each claim, margin and verdict rule was fixed in a test register before the run, and the outcome is reported whether or not it held. The quality-filter prediction held (within 0.8 points on all 12; FIGARO's own predicted retention was off by up to 7.7 points). In the development rounds the truncation recommendation was within 0.3 points of the best setting of a grid of DADA2 runs on five datasets. Soil V4 is the exception. The recommendation was 6.8 points below the best cell with an earlier build and 4.7 points below it with the final build. In both cases the best cells give the forward read more of the total and the reverse read less. A pre-registered test of such forward-heavier splits found a gain in one dataset only, so the balanced split was kept. Merges at the wrong length, the failure that a short overlap invites, were 0.00% at every cut tested on real soil V4 reads, against a reference of full-length merges. The yield estimate ranked settings in DADA2's order on 9 of 11 datasets.

Some results were negative, and they are reported. The `maxEE` rule did not meet its pre-registered adoption criteria. On the five datasets where it relaxes `maxEE`, it merged more pairs than the default on all five (by 6.1 to 50.4 points). In the pre-registered test, with a cap of 5 that has since been removed, it merged about as many pairs as relaxing to `2,5` on average (60.7% against 61.2%), but fewer on three of the five datasets; the released rule, measured afterwards on the same datasets, merged 3.7 points more on average (64.9% against 61.2%) and fewer on two. The criteria were not applied again to these numbers. It has not been tested on independent data. On the one ITS dataset the filter rows of the report matched DADA2's `filterAndTrim` exactly (74.30% of 49,720 pairs passing), while the estimate of merging was optimistic by about 10 points, with a further 15.9 points lost in denoising. These are described in the validation archive, with the scripts to reproduce them.

AmpliTrim is at version 0.1.0 and has no external users to report. [AUTHOR TO CONFIRM OR EXTEND: if AmpliTrim is already used in your own pipelines or by collaborators, say so here, with specifics.]

# AI usage disclosure

[AUTHOR TO CONFIRM AND EDIT: this paragraph must be true and complete, and must be yours.] The code, the tests, the documentation and the first draft of this paper were produced in extended working sessions with Claude, an AI assistant from Anthropic. The author set the research questions and the design. The author agreed the criteria of each validation before it was run: the assistant drafted them in a test register, and the author read and approved them before the run. The author ran all analyses on their own machine with the public data, and made every decision about what the software recommends, including keeping the `maxEE` rule after it did not meet its criteria. The author reviewed and edited the software and this paper, and takes responsibility for them.

# Acknowledgements

[AUTHOR TO ADD: funding and thanks, or delete this section.] The reads used for validation were produced by the studies whose public accessions are listed in the validation archive.

# References
