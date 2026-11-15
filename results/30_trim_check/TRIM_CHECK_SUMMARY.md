# AmpliTrim 0.3: the trim check (summary)

AmpliTrim 0.3.0; the criteria (P115, P116) are in TEST_REGISTER.md and were fixed before this ran.

| Dataset | Pairs analysed | Pass at 2,2: predicted / DADA2 | Pass after the trim: predicted / DADA2 | Yield after the trim: predicted / DADA2 merged |
|---|---|---|---|---|
| atcc_its1 | 14,898 | 74.34% / 74.34% | 91.51% / 91.51% | 86.79% / 90.50% |
| atcc_its2 | 14,818 | 26.93% / 26.93% | 91.28% / 91.28% | 80.77% / 83.33% |
| arfken_its1 | 37,076 | 61.99% / 62.20% | 90.57% / 90.66% | 83.95% / 71.06% |
| arfken_its2 | 36,721 | 37.61% / 37.70% | 90.83% / 90.91% | 85.02% / 89.15% |
| arfken_mix_its1 | 44,547 | 74.49% / 74.49% | 92.38% / 92.38% | 86.80% / 88.92% |
| arfken_mix_its2 | 44,269 | 38.06% / 38.06% | 92.48% / 92.48% | 86.20% / 86.27% |
| its2 | 43,405 | 82.86% / 82.86% | 86.87% / 86.87% | 84.50% / 14.35% |
| its | 44,665 | 74.32% / 74.24% | 94.55% / 94.44% | 87.05% / 59.54% |
| bakker_pcr_its1 | 42,294 | 91.14% / 91.04% | 96.68% / 96.66% | 91.91% / 92.35% |
| bakker_pol_its1 | 41,728 | 93.53% / 93.53% | 98.84% / 98.82% | 97.85% / 96.11% |
| bakker_pol_its2 | 37,048 | 91.46% / 91.39% | 99.06% / 98.98% | 95.50% / 96.41% |

- **P115** MET: the predicted pass rate after the trim is within 1 point of DADA2's on every dataset: 11 datasets, the largest +0.11 points (its)
- **P116** MET: the predicted pass rate at maxEE 2,2 is within 1 point of DADA2's on every dataset: 11 datasets, the largest -0.21 points (arfken_its1)

Next: send this file.
