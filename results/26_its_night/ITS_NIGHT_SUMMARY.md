# ITS overnight run: summary

The criteria (P83 to P93) are in TEST_REGISTER.md and were fixed before this ran.

- **P83** MET: On the second ITS study, at least one stricter estimate is within 10 points of DADA2's merged share: 10 differences +68.17, 5 differences +66.18, 2 differences +48.26, 0 differences -2.75 (closest: 0)
- **P84** NOT MET: On the 12 16S controls, the estimate with 0 differences is not farther from DADA2 than today's (10) by more than 2 points: mean absolute difference 8.15 against 5.84 points
- **P85 ee23** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +5.6 points; ITS2 spurious +0.00 points, recall -0.0 points
- **P86 ee23** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +9.6 points
- **P87 ee23** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +4%, its +0%
- **P85 ee24** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +11.1 points; ITS2 spurious +0.00 points, recall -0.0 points
- **P86 ee24** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +13.6 points
- **P87 ee24** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +4%, its +0%
- **P85 ee25** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +11.1 points; ITS2 spurious +0.00 points, recall -0.0 points
- **P86 ee25** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +18.0 points
- **P87 ee25** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +4%, its +0%
- **P85 ee33** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +5.6 points; ITS2 spurious +0.00 points, recall -0.0 points
- **P86 ee33** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +10.2 points
- **P87 ee33** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +9%, its +0%
- **P85 ee55** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +11.1 points; ITS2 spurious +0.00 points, recall -0.0 points
- **P86 ee55** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +20.7 points
- **P87 ee55** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +9%, its +1%
- **P85 tq10** NOT MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall -0.0 points; ITS2 spurious +0.00 points, recall -11.1 points
- **P86 tq10** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +50.5 points
- **P87 tq10** NOT MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +39%, its +0%
- **P85 tq20** NOT MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall -11.1 points; ITS2 spurious +2.38 points, recall -5.6 points
- **P86 tq20** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +20.8 points
- **P87 tq20** NOT MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +39%, its -67%
- **P85 sw15** MET: safe against std on the ATCC standard: ITS1 spurious +0.00 points, recall +16.7 points; ITS2 spurious +0.00 points, recall +27.8 points
- **P86 sw15** MET: on the ATCC standard's ITS2 the yield rises by at least 5 points: +57.3 points
- **P87 sw15** MET: a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: its2 +4%, its -2%
- **P88** MET: forward reads alone: spurious at most 1 point above std and recall not below it, on the ATCC standard: ITS1 spurious +0.15 points, recall +11.1 points; ITS2 spurious +0.57 points, recall +0.0 points
- **P89** MET: merging first: spurious at most 1 point above std and recall not below it, on the ATCC standard: ITS1 spurious +0.00 points, recall +16.7 points; ITS2 spurious +0.24 points, recall +33.3 points
- **P90** MET: merging first raises the yield of the second ITS study by at least 10 points: +50.4 points
- **P91** MET: AmpliTrim's predicted yields order the maxEE arms as DADA2's merged shares do, rank correlation at least 0.8 on each ITS group: atcc_its1 1.00 (6 arms), atcc_its2 1.00 (6 arms), arfken_its1 1.00 (6 arms), arfken_its2 1.00 (6 arms), arfken_mix_its1 1.00 (6 arms), arfken_mix_its2 1.00 (6 arms), its2 1.00 (6 arms), its 1.00 (6 arms)
- **P92** MET: AmpliTrim's predicted pass rate is within 1 point of DADA2's at every maxEE arm on every ITS group: 48 comparisons, the largest +0.22 points (its, ee25)
- **P93** MET: the spurious share at 12 nt is at most 0.5 point above 20 and 30 nt on the ATCC standard: ITS1 12 nt against 20 nt +0.00, against 30 nt +0.00 points; ITS2 12 nt against 20 nt +0.00, against 30 nt +0.00 points

**The advice for ITS (section 33):** `sw15`, the arm that is safe on both markers and useful on ITS2 with the largest ITS2 yield.

Next: send this file.
