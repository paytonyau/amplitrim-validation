# Recovery of the known members of the mock communities

Members found per sample: holding at least 0.1% of the sample's reads (the rule of sections 33 to 44), and holding at least one read. Then the members found in any sample and in all samples (of those expected in every sample). ITS: run 31 (std: DADA2's defaults; sw15: the trim). 16S: run 26 (DADA2 at AmpliTrim's cut, maxEE 2,2), a strain found when one of its variants is (exact, as mock_check.py).

| Marker | Dataset | Arm | Samples | Expected per sample | Found per sample, at least 0.1% | Found per sample, any read | In any sample (0.1% / any read) | In all samples (0.1% / any read) |
|---|---|---|---|---|---|---|---|---|
| ITS | atcc_its1 | std | 3 | 6.0 | 5.0 (83.3%; 5 to 5) | 5.3 (88.9%) | 5 / 6 of 6 | 5 / 5 of 6 |
| ITS | atcc_its1 | sw15 | 3 | 6.0 | 6.0 (100.0%; 6 to 6) | 6.0 (100.0%) | 6 / 6 of 6 | 6 / 6 of 6 |
| ITS | atcc_its2 | std | 3 | 6.0 | 4.0 (66.7%; 4 to 4) | 4.0 (66.7%) | 4 / 4 of 6 | 4 / 4 of 6 |
| ITS | atcc_its2 | sw15 | 3 | 6.0 | 5.7 (94.4%; 5 to 6) | 5.7 (94.4%) | 6 / 6 of 6 | 5 / 5 of 6 |
| ITS | arfken_mix_its1 | std | 9 | 3.0 | 2.3 (76.7%; 1 to 4) | 2.3 (76.7%) | 4 / 4 of 5 | 1 / 1 of 1 |
| ITS | arfken_mix_its1 | sw15 | 9 | 3.0 | 2.3 (76.7%; 1 to 4) | 2.3 (76.7%) | 4 / 4 of 5 | 1 / 1 of 1 |
| ITS | arfken_mix_its2 | std | 9 | 3.0 | 3.0 (100.0%; 2 to 5) | 3.0 (100.0%) | 5 / 5 of 5 | 1 / 1 of 1 |
| ITS | arfken_mix_its2 | sw15 | 9 | 3.0 | 3.0 (100.0%; 2 to 5) | 3.0 (100.0%) | 5 / 5 of 5 | 1 / 1 of 1 |
| ITS | bakker_pcr_its1 | std | 27 | 18.0 | 10.6 (58.8%; 5 to 14) | 11.1 (61.5%) | 14 / 14 of 18 | 5 / 6 of 18 |
| ITS | bakker_pcr_its1 | sw15 | 27 | 18.0 | 10.5 (58.2%; 6 to 14) | 10.7 (59.7%) | 14 / 14 of 18 | 5 / 6 of 18 |
| ITS | bakker_pol_its1 | std | 27 | 18.0 | 9.4 (52.1%; 5 to 14) | 9.7 (53.7%) | 17 / 17 of 18 | 5 / 5 of 18 |
| ITS | bakker_pol_its1 | sw15 | 27 | 18.0 | 9.8 (54.3%; 5 to 16) | 10.1 (56.2%) | 18 / 18 of 18 | 5 / 5 of 18 |
| ITS | bakker_pol_its2 | std | 27 | 19.0 | 11.7 (61.4%; 7 to 15) | 11.7 (61.8%) | 15 / 15 of 19 | 5 / 5 of 19 |
| ITS | bakker_pol_its2 | sw15 | 27 | 19.0 | 11.4 (60.2%; 7 to 15) | 11.5 (60.6%) | 15 / 15 of 19 | 5 / 5 of 19 |
| ITS | pauvert_its1 | std | 3 | 151.0 | 87.3 (57.8%; 86 to 89) | 122.7 (81.2%) | 97 / 133 of 151 (155 strains) | 81 / 114 of 151 |
| ITS | pauvert_its1 | sw15 | 3 | 151.0 | 85.7 (56.7%; 84 to 87) | 122.0 (80.8%) | 93 / 133 of 151 (155 strains) | 80 / 112 of 151 |
| 16S | mockA_high | std | 6 | 14.0 | 12.0 (85.7%; 12 to 12) | 12.0 (85.7%) | 12 / 12 of 14 | 12 / 12 of 14 |
| 16S | mockB | std | 8 | 3.0 | 2.9 (95.8%; 2 to 3) | 2.9 (95.8%) | 3 / 3 of 3 | 2 / 2 of 3 |
| 16S | d6300_V3-V4 | std | 10 | 8.0 | 8.0 (100.0%; 8 to 8) | 8.0 (100.0%) | 8 / 8 of 8 | 8 / 8 of 8 |
| 16S | d6300_V4-V5 | std | 3 | 9.0 | 8.0 (88.9%; 8 to 8) | 8.0 (88.9%) | 8 / 8 of 9 | 8 / 8 of 9 |
| 16S | miseqsop_mock | std | 1 | 21.0 | 19.0 (90.5%; 19 to 19) | 19.0 (90.5%) | 19 / 19 of 21 | 19 / 19 of 21 |

Published, for comparison: pauvert_its1: 75.8% of the 189 strains per replicate (Pauvert et al. 2019, their recommended DADA2 approach, Si5).
16S exact variants recovered (pooled): mockA_high 36 of 53; mockB 6 of 7; d6300_V3-V4 10 of 10; d6300_V4-V5 9 of 10; miseqsop_mock 20 of 32.

Next: send this file.
