# ITS rules night: summary

The criteria (P94 to P104) are in TEST_REGISTER.md and were fixed before this ran.

**The consistency check:** std reproduces run 26 on 8 datasets exactly.

- **P94 tq3** MET: safe on every truth: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P95 tq3** NOT MET: the ATCC standard's ITS2 yield rises by at least 5 points: +0.0
- **P96 tq3** NOT MET: the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: +0.00
- **P97 tq3** MET: a flag: the median ASVs per sample rise by at most 25%: its +0%
- **P94 tq4** MET: safe on every truth: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P95 tq4** NOT MET: the ATCC standard's ITS2 yield rises by at least 5 points: +0.0
- **P96 tq4** NOT MET: the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: +0.00
- **P97 tq4** MET: a flag: the median ASVs per sample rise by at most 25%: its +0%
- **P94 tq5** MET: safe on every truth: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P95 tq5** NOT MET: the ATCC standard's ITS2 yield rises by at least 5 points: +0.0
- **P96 tq5** NOT MET: the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: +0.00
- **P97 tq5** MET: a flag: the median ASVs per sample rise by at most 25%: its +0%
- **P94 tq6** MET: safe on every truth: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P95 tq6** NOT MET: the ATCC standard's ITS2 yield rises by at least 5 points: +0.0
- **P96 tq6** NOT MET: the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: +0.00
- **P97 tq6** MET: a flag: the median ASVs per sample rise by at most 25%: its +0%
- **P94 tq8** NOT MET: safe on every truth: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -2.2
- **P95 tq8** MET: the ATCC standard's ITS2 yield rises by at least 5 points: +54.5
- **P96 tq8** NOT MET: the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: -1.08
- **P97 tq8** MET: a flag: the median ASVs per sample rise by at most 25%: its +0%
- **P98 sw10** MET: safe on every truth: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +5.6; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P99 sw10** MET: the ATCC standard's ITS2 yield rises by at least 5 points: +18.7
- **P100 sw10** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +4%, its +0%
- **P98 sw20** NOT MET: safe on every truth: atcc_its1 spurious +0.00, recall +11.1; atcc_its2 spurious +0.00, recall +11.1; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -4.4
- **P99 sw20** MET: the ATCC standard's ITS2 yield rises by at least 5 points: +57.7
- **P100 sw20** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its -8%
- **P98 sw25** NOT MET: safe on every truth: atcc_its1 spurious +0.00, recall -0.0; atcc_its2 spurious +5.99, recall -0.0; arfken_mix_its1 spurious -0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -4.4
- **P99 sw25** MET: the ATCC standard's ITS2 yield rises by at least 5 points: +38.7
- **P100 sw25** MET: a flag: the median ASVs per sample rise by at most 25%: its2 -4%, its -28%
- **P101** MET: sw15 reproduces run 26's sw15 on 8 datasets
- **P102** MET: safe on every truth: atcc_its1 spurious +0.00, recall +5.6; atcc_its2 spurious +0.00, recall -0.0; arfken_mix_its1 spurious +0.00, recall -0.0; arfken_mix_its2 spurious +0.00, recall -0.0
- **P103** MET: the ATCC standard's ITS2 yield rises by at least 5 points: +33.6
- **P104** MET: a flag: the median ASVs per sample rise by at most 25%: its2 +0%, its +0%

**What follows (section 35):** no truncQ value from 3 to 8 is safe and useful: truncQ stays 2. The trim strength: `sw15`. Advised for poor read ends: `sw15`. The computed maxEE (`eer`) may be used for ITS.

Next: send this file.
