# The v0.4 night: summary

The criteria (P117 to P130) are in TEST_REGISTER.md and were fixed before this ran.

**Gate, consistency** with run 29: passed (std: 11 datasets compared, identical; sw15: 11 datasets compared, identical; eer: 11 datasets compared, identical).
**Gate, overhangs** of at least 3 bases in std: bakker_pcr_its1 0.00%, bakker_pol_its1 0.00%, bakker_pol_its2 0.00%, pauvert_its1 0.00% (at most 1%: passed).
**Gate, Pauvert's truth:** 151 expected members (at least 120: passed).

- **P117** MET: pauvert_its1 spurious -0.02, recall -1.1
- **P118** MET: pauvert_its1 spurious +0.00, recall +0.0
- **P119 w2** NOT MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.03, recall -3.1; bakker_pol_its1 spurious +0.00, recall +4.5; bakker_pol_its2 spurious +0.00, recall -11.9; pauvert_its1 spurious -0.21, recall -1.5
- **P120 w2** NOT MET: ATCC ITS2 yield -0.0 points against sw15 (at least +2), lowest -18.4 (at least -1)
- **P119 w3** NOT MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall -2.2; bakker_pcr_its1 spurious -0.03, recall -1.4; bakker_pol_its1 spurious +0.00, recall +2.5; bakker_pol_its2 spurious +0.00, recall -1.9; pauvert_its1 spurious -0.03, recall -1.3
- **P120 w3** NOT MET: ATCC ITS2 yield +0.6 points against sw15 (at least +2), lowest -0.6 (at least -1)
- **P119 w5** MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +27.8; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.03, recall -0.4; bakker_pol_its1 spurious +0.00, recall +4.1; bakker_pol_its2 spurious +0.00, recall -0.2; pauvert_its1 spurious +0.00, recall -0.4
- **P120 w5** NOT MET: ATCC ITS2 yield -3.9 points against sw15 (at least +2), lowest -3.9 (at least -1)
- **P119 w6** NOT MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +16.7; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall -2.2; bakker_pcr_its1 spurious -0.03, recall +0.6; bakker_pol_its1 spurious +0.00, recall +3.9; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious +0.03, recall -0.2
- **P120 w6** NOT MET: ATCC ITS2 yield -10.4 points against sw15 (at least +2), lowest -10.4 (at least -1)
- **P119 w8** MET: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +22.2; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.02, recall +0.6; bakker_pol_its1 spurious +0.00, recall +1.9; bakker_pol_its2 spurious +0.00, recall -0.4; pauvert_its1 spurious +0.04, recall +0.0
- **P120 w8** NOT MET: ATCC ITS2 yield -19.0 points against sw15 (at least +2), lowest -19.0 (at least -1)
- **P119 w10** MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +22.2; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.02, recall +0.6; bakker_pol_its1 spurious +0.00, recall +1.9; bakker_pol_its2 spurious +0.00, recall -0.2; pauvert_its1 spurious +0.05, recall -0.4
- **P120 w10** NOT MET: ATCC ITS2 yield -26.7 points against sw15 (at least +2), lowest -26.7 (at least -1)
- **P119 w12** MET: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.01, recall +0.2; bakker_pol_its1 spurious +0.00, recall +1.9; bakker_pol_its2 spurious +0.00, recall -0.2; pauvert_its1 spurious +0.04, recall +0.2
- **P120 w12** NOT MET: ATCC ITS2 yield -32.7 points against sw15 (at least +2), lowest -32.7 (at least -1)
- **P119 w16** MET: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +16.7; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.01, recall +0.0; bakker_pol_its1 spurious +0.00, recall +1.4; bakker_pol_its2 spurious +0.00, recall -0.2; pauvert_its1 spurious +0.03, recall +0.2
- **P120 w16** NOT MET: ATCC ITS2 yield -37.9 points against sw15 (at least +2), lowest -37.9 (at least -1)
- **P119 w20** MET: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.01, recall +0.2; bakker_pol_its1 spurious +0.00, recall +1.0; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious +0.01, recall +0.0
- **P120 w20** NOT MET: ATCC ITS2 yield -42.2 points against sw15 (at least +2), lowest -42.2 (at least -1)
- **P119 q12** MET: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.02, recall +0.2; bakker_pol_its1 spurious +0.00, recall +0.0; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious +0.03, recall -0.2
- **P120 q12** NOT MET: ATCC ITS2 yield -20.2 points against sw15 (at least +2), lowest -20.2 (at least -1)
- **P119 q18** NOT MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall -2.2; bakker_pcr_its1 spurious -0.04, recall -2.5; bakker_pol_its1 spurious +0.00, recall +4.5; bakker_pol_its2 spurious +0.00, recall -5.5; pauvert_its1 spurious -0.05, recall -1.1
- **P120 q18** NOT MET: ATCC ITS2 yield +0.6 points against sw15 (at least +2), lowest -3.0 (at least -1)
- **P121 ml20** MET: atcc_its1 spurious +0.00, recall +0.0; atcc_its2 spurious +0.00, recall +0.0; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious +0.00, recall +0.0; bakker_pol_its1 spurious +0.00, recall +0.0; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious -0.01, recall -0.2
- **P121 ml100** MET: atcc_its1 spurious +0.00, recall +0.0; atcc_its2 spurious +0.00, recall +0.0; arfken_mix_its1 spurious +0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.00, recall -0.4; bakker_pol_its1 spurious +0.00, recall +0.0; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious +0.02, recall +0.4
- **P121 ml150** NOT MET: atcc_its1 spurious +0.00, recall +0.0; atcc_its2 spurious +0.00, recall +0.0; arfken_mix_its1 spurious +0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.96, recall -8.6; bakker_pol_its1 spurious +0.00, recall -1.2; bakker_pol_its2 spurious +0.00, recall +0.0; pauvert_its1 spurious -0.02, recall +1.1
- **P122 tq3** not decided: changes nothing on any truth dataset (not judged)
- **P123 tq3** not decided: changes nothing on any truth dataset (not judged)
- **P122 tq5** not decided: changes nothing on any truth dataset (not judged)
- **P123 tq5** not decided: changes nothing on any truth dataset (not judged)
- **P122 tq7** MET: atcc_its1 spurious +0.00, recall +16.7; atcc_its2 spurious +0.00, recall +33.3; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.01, recall -0.8; pauvert_its1 spurious -0.09, recall +0.0
- **P123 tq7** MET: the largest yield gain where it differs: +45.9 points (atcc_its1, atcc_its2, arfken_mix_its1, arfken_mix_its2, bakker_pcr_its1, pauvert_its1)
- **P124** MET: atcc_its1 spurious +0.15, recall +11.1; atcc_its2 spurious +0.57, recall +0.0; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.00, recall +5.8; bakker_pol_its1 spurious +0.03, recall +10.9; bakker_pol_its2 spurious +0.00, recall +7.4; pauvert_its1 spurious -0.99, recall +1.8; recall not below std's everywhere: True
- **P129** NOT MET: atcc_its1 spurious +3.60, recall +5.6; atcc_its2 spurious +7.22, recall +27.8; arfken_mix_its1 spurious -0.00, recall +0.0; arfken_mix_its2 spurious +0.00, recall +0.0; bakker_pcr_its1 spurious -0.03, recall +1.9; bakker_pol_its1 spurious -0.00, recall +4.7; bakker_pol_its2 spurious +0.00, recall -0.2; pauvert_its1 spurious +0.11, recall -0.7; ASVs per member found against fwd: atcc_its1 +5%, atcc_its2 +18%, arfken_mix_its1 +0%, arfken_mix_its2 +0%, bakker_pcr_its1 +5%, bakker_pol_its1 +1%, bakker_pol_its2 +0%, pauvert_its1 -1% (at most +25%)
- **P125 std** NOT MET: mean absolute difference from DADA2's merged share: 0 differences 18.1 points, 10 differences 9.4; largest at 0: 43.8 (at most 10)
- **P125 sw15** NOT MET: mean absolute difference from DADA2's merged share: 0 differences 28.1 points, 10 differences 10.6; largest at 0: 43.1 (at most 10)
- **P126** MET: the largest difference +0.25 points (its, window 10)
- **P127** MET: 12 16S datasets: d6300_V3-V4: same (224,218, maxEE 2.0,2.0); d6300_V4-V5: same (282,280, maxEE 2.0,3.0); mcm_v12: same (160,165, maxEE 2.0,2.0); mcm_v46: same (273,265, maxEE 4.0,3.0); miseqsop_mock: same (132,134, maxEE 2.0,2.0); mockA_high: same (273,247, maxEE 2.0,3.0); mockB: same (274,253, maxEE 4.0,8.0); spikein_V3-V4: same (229,210, maxEE 2.0,2.0); swabs_V3-V4: same (224,218, maxEE 2.0,2.0); v1v2: same (195,168, maxEE 2.0,2.0); v4: same (124,115, maxEE 2.0,2.0); wastewater_V3-V4: same (260,182, maxEE 4.0,3.0)
- **P128 sw15** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -2%
- **P128 w2** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its -12%
- **P128 w3** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -5%
- **P128 w5** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its +1%
- **P128 w6** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +1%
- **P128 w8** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +1%
- **P128 w10** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +2%
- **P128 w12** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +2%
- **P128 w16** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +2%
- **P128 w20** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +2%
- **P128 q12** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its +0%
- **P128 q18** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its -4%
- **P128 eer** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P128 tq3** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P128 tq5** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P128 tq7** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%
- **P128 fwd** NOT MET: a flag: the median ASVs per sample rise by at most 25%: its2 +309%, its +30%
- **P130** MET: the best trim: sw15 (tied within 0.5 points: sw15)

## The trims (B)

| Trim | Safe on every truth | Mean yield (truth) | ATCC ITS2 yield | Mean bases cut per read | Mean merge share |
|---|---|---|---|---|---|
| sw15 (best) | True | 88.53% | 83.01% | 23.9 | 96.10% |
| w2 | False | 84.75% | 82.99% | 40.9 | 94.34% |
| w3 | False | 88.64% | 83.66% | 30.2 | 95.78% |
| w5 | True | 87.02% | 79.07% | 20.1 | 96.37% |
| w6 | False | 84.61% | 72.61% | 17.7 | 96.28% |
| w8 | True | 80.61% | 63.98% | 14.2 | 96.66% |
| w10 | True | 77.72% | 56.27% | 12.0 | 96.94% |
| w12 | True | 75.69% | 50.33% | 10.5 | 96.92% |
| w16 | True | 73.53% | 45.14% | 8.8 | 97.54% |
| w20 | True | 72.10% | 40.81% | 7.8 | 97.80% |
| q12 | True | 80.09% | 62.79% | 14.5 | 96.45% |
| q18 | False | 87.59% | 83.57% | 33.8 | 95.76% |

Next: send this file.
