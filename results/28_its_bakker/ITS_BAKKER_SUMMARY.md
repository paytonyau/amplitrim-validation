# The Bakker night: summary

The criteria (P105 to P114) are in TEST_REGISTER.md and were fixed before this ran.

**The consistency check:** std reproduces run 26 on 8 datasets exactly.

- **P105** MET: the trim (sw15) is safe on every truth: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious -4.41, recall -1.0; bakker_pol_its1 spurious +0.85, recall +2.1; bakker_pol_its2 spurious +0.00, recall -1.2
- **P106** MET: on every truth dataset where fewer than 50% of the pairs pass at maxEE 2,2, the trim raises the yield by at least 5 points: atcc_its2 +57.3, arfken_mix_its2 +46.7
- **P107** MET: the 50% trigger: where at least 50% of the pairs pass at maxEE 2,2, the trim is safe and loses at most 2 points of yield: atcc_its1 safe, yield +17.0; arfken_mix_its1 safe, yield +15.4; bakker_pcr_its1 safe, yield -1.4; bakker_pol_its1 safe, yield +5.3; bakker_pol_its2 safe, yield +6.2
- **P108** MET: the computed maxEE (eer) is safe on every truth: atcc_its1 spurious +0.00, recall +5.6; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious +0.00, recall -0.0; bakker_pol_its1 spurious +0.00, recall -0.0; bakker_pol_its2 spurious +0.00, recall -0.0
- **P109** MET: the computed maxEE loses at most 1 point of yield on every truth dataset: atcc_its1 +4.5, atcc_its2 +33.6, arfken_mix_its1 +0.0, arfken_mix_its2 +29.0, bakker_pcr_its1 +0.0, bakker_pol_its1 +0.0, bakker_pol_its2 +0.0
- **P110** MET: the trim with maxEE computed on the trimmed reads (sweer) is safe on every truth: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0; bakker_pcr_its1 spurious -4.41, recall -1.0; bakker_pol_its1 spurious +0.85, recall +2.1; bakker_pol_its2 spurious +0.00, recall -1.2
- **P111** NOT MET: sweer against the trim alone: the ATCC standard's ITS2 yield +0.0 points (at least +2), and nowhere lower than -1: atcc_its1 +0.0, atcc_its2 +0.0, arfken_mix_its1 +0.0, arfken_mix_its2 +0.0, bakker_pcr_its1 +0.0, bakker_pol_its1 +0.0, bakker_pol_its2 +0.0
- **P112 sw15** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -2%
- **P112 eer** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P112 sweer** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -2%
- **P113** MET: sw15 reproduces run 26 on 8 datasets and eer reproduces run 27 on 8
- **P114** MET: the known-unsafe reference (truncQ 10) is still unsafe somewhere: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -11.1; arfken_mix_its1 spurious -0.00, recall -30.0; arfken_mix_its2 spurious +0.00, recall -2.2; bakker_pcr_its1 spurious -1.57, recall -1.0; bakker_pol_its1 spurious +0.00, recall -0.0; bakker_pol_its2 spurious +0.00, recall -0.0

**What follows (section 37):** The trim advice (window 4, mean quality 15) holds on a third mock community, with other read lengths and error-prone PCR conditions. The trim did not harm data that pass well at maxEE 2,2: the 50% trigger is not needed for safety; it may stay as a guide to when the trim helps. The computed maxEE for ITS holds on a third mock community. The combination is safe but not better than the trim alone: the trim and the computed maxEE stay separate.

Next: send this file.
