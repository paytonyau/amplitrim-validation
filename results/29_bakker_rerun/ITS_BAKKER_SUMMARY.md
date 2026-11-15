# The Bakker night: summary

The criteria (P105 to P114) are in TEST_REGISTER.md and were fixed before this ran.

**The consistency check:** std reproduces run 26 on 8 datasets exactly.

**The overhang check:** reads of std carrying primer overhangs: bakker_pcr_its1 1.69%, bakker_pol_its1 31.40%, bakker_pol_its2 0.00%: FAILED (more than 1%, or a dataset missing), so no criterion is decided.

- **P105** not decided (the consistency check failed): the trim (sw15) is safe on every truth: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious -0.04, recall -0.6; bakker_pol_its1 spurious +0.00, recall +2.3; bakker_pol_its2 spurious +0.00, recall -1.2
- **P106** not decided (the consistency check failed): on every truth dataset where fewer than 50% of the pairs pass at maxEE 2,2, the trim raises the yield by at least 5 points: atcc_its2 +57.3, arfken_mix_its2 +46.7
- **P107** not decided (the consistency check failed): the 50% trigger: where at least 50% of the pairs pass at maxEE 2,2, the trim is safe and loses at most 2 points of yield: atcc_its1 safe, yield +17.0; arfken_mix_its1 safe, yield +15.4; bakker_pcr_its1 safe, yield +5.2; bakker_pol_its1 safe, yield +5.1; bakker_pol_its2 safe, yield +6.2
- **P108** not decided (the consistency check failed): the computed maxEE (eer) is safe on every truth: atcc_its1 spurious +0.00, recall +5.6; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious +0.00, recall -0.0; bakker_pol_its1 spurious +0.00, recall -0.0; bakker_pol_its2 spurious +0.00, recall -0.0
- **P109** not decided (the consistency check failed): the computed maxEE loses at most 1 point of yield on every truth dataset: atcc_its1 +4.5, atcc_its2 +33.6, arfken_mix_its1 +0.0, arfken_mix_its2 +29.0, bakker_pcr_its1 +0.0, bakker_pol_its1 +0.0, bakker_pol_its2 +0.0
- **P110** not decided (the consistency check failed): the trim with maxEE computed on the trimmed reads (sweer) is safe on every truth: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious -0.04, recall -0.6; bakker_pol_its1 spurious +0.00, recall +2.3; bakker_pol_its2 spurious +0.00, recall -1.2
- **P111** not decided (the consistency check failed): sweer against the trim alone: the ATCC standard's ITS2 yield +0.0 points (at least +2), and nowhere lower than -1: atcc_its1 +0.0, atcc_its2 +0.0, arfken_mix_its1 +0.0, arfken_mix_its2 +0.0, bakker_pcr_its1 +0.0, bakker_pol_its1 +0.0, bakker_pol_its2 +0.0
- **P112 sw15** not decided (the consistency check failed): a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -2%
- **P112 eer** not decided (the consistency check failed): a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P112 sweer** not decided (the consistency check failed): a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -2%
- **P113** not decided (the consistency check failed): sw15 reproduces run 26 on 8 datasets and eer reproduces run 27 on 8
- **P114** not decided (the consistency check failed): the known-unsafe reference (truncQ 10) is still unsafe somewhere: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -11.1; arfken_mix_its1 spurious -0.00, recall -30.0; arfken_mix_its2 spurious +0.00, recall -2.2; bakker_pcr_its1 spurious -0.01, recall -4.1; bakker_pol_its1 spurious +0.00, recall -0.0; bakker_pol_its2 spurious +0.00, recall -0.0

**What follows (section 37):** No criterion was decided.

Next: send this file.
