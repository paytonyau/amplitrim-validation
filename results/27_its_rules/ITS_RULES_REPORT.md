# ITS rules night: report

## The computed maxEE per dataset (the arm eer)

| Dataset | P75 of EE, forward / reverse | maxEE | Reads used |
|---|---|---|---|
| arfken_its1 | 1.02 / 3.94 | 2,4 | 74,151 / 74,151 |
| arfken_its2 | 1.85 / 5.74 | 2,6 | 73,441 / 73,441 |
| arfken_mix_its1 | 0.57 / 1.95 | 2,2 | 44,547 / 44,547 |
| arfken_mix_its2 | 1.82 / 5.53 | 2,6 | 44,269 / 44,269 |
| atcc_its1 | 0.70 / 2.02 | 2,3 | 14,898 / 14,898 |
| atcc_its2 | 2.99 / 8.00 | 3,9 | 14,818 / 14,818 |
| its | 0.89 / 1.73 | 2,2 | 178,657 / 178,657 |
| its2 | 0.84 / 0.80 | 2,2 | 43,405 / 43,405 |

## Flows and merging (shares of the input pairs; merge share: of the denoised pairs)

| Dataset | Arm | Input | Pass | Merged | After chimera removal | Merge share | ASVs |
|---|---|---|---|---|---|---|---|
| arfken_its1 | eer | 74,151 | 73.02% | 59.34% | 59.07% | 81.39% | 2 |
| arfken_its1 | std | 74,151 | 62.20% | 52.60% | 52.52% | 84.69% | 2 |
| arfken_its1 | sw10 | 74,151 | 70.50% | 58.39% | 58.12% | 82.95% | 2 |
| arfken_its1 | sw15 | 74,151 | 90.66% | 71.06% | 70.59% | 78.67% | 2 |
| arfken_its1 | sw20 | 74,151 | 85.47% | 66.62% | 66.23% | 78.26% | 2 |
| arfken_its1 | sw25 | 74,151 | 77.92% | 60.55% | 60.46% | 77.90% | 2 |
| arfken_its1 | tq10 | 74,151 | 91.42% | 78.02% | 77.91% | 85.72% | 1 |
| arfken_its1 | tq3 | 74,151 | 62.20% | 52.60% | 52.52% | 84.69% | 2 |
| arfken_its1 | tq4 | 74,151 | 62.20% | 52.60% | 52.52% | 84.69% | 2 |
| arfken_its1 | tq5 | 74,151 | 62.20% | 52.60% | 52.52% | 84.69% | 2 |
| arfken_its1 | tq6 | 74,151 | 62.20% | 52.60% | 52.52% | 84.69% | 2 |
| arfken_its1 | tq8 | 74,151 | 84.40% | 82.00% | 81.59% | 97.35% | 2 |
| arfken_its2 | eer | 73,441 | 70.26% | 69.89% | 69.89% | 99.66% | 1 |
| arfken_its2 | std | 73,441 | 37.70% | 37.57% | 37.57% | 99.73% | 1 |
| arfken_its2 | sw10 | 73,441 | 56.40% | 56.16% | 56.16% | 99.67% | 1 |
| arfken_its2 | sw15 | 73,441 | 90.91% | 89.15% | 89.15% | 98.32% | 1 |
| arfken_its2 | sw20 | 73,441 | 88.96% | 85.68% | 85.68% | 96.49% | 1 |
| arfken_its2 | sw25 | 73,441 | 79.20% | 76.60% | 76.60% | 96.89% | 1 |
| arfken_its2 | tq10 | 73,441 | 90.01% | 86.90% | 86.90% | 96.83% | 1 |
| arfken_its2 | tq3 | 73,441 | 37.70% | 37.57% | 37.57% | 99.73% | 1 |
| arfken_its2 | tq4 | 73,441 | 37.70% | 37.57% | 37.57% | 99.73% | 1 |
| arfken_its2 | tq5 | 73,441 | 37.70% | 37.57% | 37.57% | 99.73% | 1 |
| arfken_its2 | tq6 | 73,441 | 37.70% | 37.57% | 37.57% | 99.73% | 1 |
| arfken_its2 | tq8 | 73,441 | 86.46% | 84.23% | 84.23% | 97.61% | 1 |
| arfken_mix_its1 | eer | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | std | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | sw10 | 44,547 | 78.95% | 76.75% | 75.33% | 97.33% | 5 |
| arfken_mix_its1 | sw15 | 44,547 | 92.38% | 88.92% | 87.52% | 96.52% | 4 |
| arfken_mix_its1 | sw20 | 44,547 | 86.98% | 84.04% | 82.76% | 96.79% | 4 |
| arfken_mix_its1 | sw25 | 44,547 | 79.71% | 77.35% | 76.57% | 97.31% | 4 |
| arfken_mix_its1 | tq10 | 44,547 | 91.66% | 80.42% | 80.14% | 88.06% | 2 |
| arfken_mix_its1 | tq3 | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | tq4 | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | tq5 | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | tq6 | 44,547 | 74.49% | 72.42% | 72.08% | 97.42% | 4 |
| arfken_mix_its1 | tq8 | 44,547 | 86.45% | 83.71% | 82.29% | 97.19% | 4 |
| arfken_mix_its2 | eer | 44,269 | 71.34% | 69.34% | 63.01% | 97.77% | 3 |
| arfken_mix_its2 | std | 44,269 | 38.06% | 37.28% | 33.97% | 98.72% | 3 |
| arfken_mix_its2 | sw10 | 44,269 | 57.67% | 55.67% | 51.21% | 97.10% | 3 |
| arfken_mix_its2 | sw15 | 44,269 | 92.48% | 86.27% | 80.67% | 93.88% | 2 |
| arfken_mix_its2 | sw20 | 44,269 | 90.49% | 86.81% | 79.20% | 96.38% | 2 |
| arfken_mix_its2 | sw25 | 44,269 | 81.58% | 80.18% | 73.91% | 98.68% | 2 |
| arfken_mix_its2 | tq10 | 44,269 | 92.29% | 88.83% | 81.14% | 96.74% | 2 |
| arfken_mix_its2 | tq3 | 44,269 | 38.06% | 37.28% | 33.97% | 98.72% | 3 |
| arfken_mix_its2 | tq4 | 44,269 | 38.06% | 37.28% | 33.97% | 98.72% | 3 |
| arfken_mix_its2 | tq5 | 44,269 | 38.06% | 37.28% | 33.97% | 98.72% | 3 |
| arfken_mix_its2 | tq6 | 44,269 | 38.06% | 37.28% | 33.97% | 98.72% | 3 |
| arfken_mix_its2 | tq8 | 44,269 | 87.72% | 84.39% | 77.05% | 96.69% | 2 |
| atcc_its1 | eer | 14,898 | 78.92% | 78.04% | 77.86% | 99.29% | 9 |
| atcc_its1 | std | 14,898 | 74.34% | 73.53% | 73.35% | 99.28% | 8 |
| atcc_its1 | sw10 | 14,898 | 77.89% | 77.10% | 76.91% | 99.32% | 10 |
| atcc_its1 | sw15 | 14,898 | 91.51% | 90.50% | 90.32% | 99.23% | 10 |
| atcc_its1 | sw20 | 14,898 | 86.80% | 85.87% | 85.70% | 99.33% | 10 |
| atcc_its1 | sw25 | 14,898 | 79.19% | 78.11% | 77.96% | 99.00% | 9 |
| atcc_its1 | tq10 | 14,898 | 91.54% | 90.05% | 89.92% | 98.89% | 9 |
| atcc_its1 | tq3 | 14,898 | 74.34% | 73.53% | 73.35% | 99.28% | 8 |
| atcc_its1 | tq4 | 14,898 | 74.34% | 73.53% | 73.35% | 99.28% | 8 |
| atcc_its1 | tq5 | 14,898 | 74.34% | 73.53% | 73.35% | 99.28% | 8 |
| atcc_its1 | tq6 | 14,898 | 74.34% | 73.53% | 73.35% | 99.28% | 8 |
| atcc_its1 | tq8 | 14,898 | 85.43% | 84.55% | 84.36% | 99.26% | 10 |
| atcc_its2 | eer | 14,818 | 70.99% | 60.02% | 59.26% | 85.98% | 7 |
| atcc_its2 | std | 14,818 | 26.93% | 26.14% | 25.70% | 98.63% | 5 |
| atcc_its2 | sw10 | 14,818 | 47.70% | 44.84% | 44.41% | 95.61% | 8 |
| atcc_its2 | sw15 | 14,818 | 91.28% | 83.33% | 83.01% | 91.95% | 10 |
| atcc_its2 | sw20 | 14,818 | 88.95% | 85.00% | 83.38% | 95.97% | 8 |
| atcc_its2 | sw25 | 14,818 | 79.11% | 65.87% | 64.37% | 83.53% | 7 |
| atcc_its2 | tq10 | 14,818 | 92.16% | 77.62% | 76.22% | 84.69% | 7 |
| atcc_its2 | tq3 | 14,818 | 26.93% | 26.14% | 25.70% | 98.63% | 5 |
| atcc_its2 | tq4 | 14,818 | 26.93% | 26.14% | 25.70% | 98.63% | 5 |
| atcc_its2 | tq5 | 14,818 | 26.93% | 26.14% | 25.70% | 98.63% | 5 |
| atcc_its2 | tq6 | 14,818 | 26.93% | 26.14% | 25.70% | 98.63% | 5 |
| atcc_its2 | tq8 | 14,818 | 84.11% | 81.63% | 80.22% | 97.55% | 9 |
| its | eer | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | std | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | sw10 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | sw15 | 178,657 | 94.44% | 59.54% | 59.51% | 80.35% | 103 |
| its | sw20 | 178,657 | 91.62% | 58.46% | 58.45% | 81.95% | 97 |
| its | sw25 | 178,657 | 83.39% | 49.53% | 49.53% | 79.39% | 76 |
| its | tq10 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | tq3 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | tq4 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | tq5 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | tq6 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its | tq8 | 178,657 | 74.24% | 46.87% | 46.85% | 80.29% | 105 |
| its2 | eer | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | std | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | sw10 | 43,405 | 83.95% | 14.30% | 14.29% | 26.77% | 24 |
| its2 | sw15 | 43,405 | 86.87% | 14.35% | 14.35% | 26.07% | 24 |
| its2 | sw20 | 43,405 | 97.00% | 15.95% | 15.92% | 25.51% | 23 |
| its2 | sw25 | 43,405 | 97.83% | 16.63% | 16.60% | 26.68% | 22 |
| its2 | tq10 | 43,405 | 91.85% | 31.69% | 31.28% | 60.17% | 32 |
| its2 | tq3 | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | tq4 | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | tq5 | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | tq6 | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |
| its2 | tq8 | 43,405 | 82.86% | 13.96% | 13.96% | 26.46% | 23 |

## Accuracy (means over samples)

| Dataset | Arm | Spurious | Recall | Exact | Members found | Unjudged |
|---|---|---|---|---|---|---|
| arfken_mix_its1 | eer | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | std | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | sw10 | 0.01% | 76.67% | 99.99% | 2.3 | 0.00% |
| arfken_mix_its1 | sw15 | 0.01% | 76.67% | 99.99% | 2.3 | 0.00% |
| arfken_mix_its1 | sw20 | 0.01% | 76.67% | 99.99% | 2.3 | 0.00% |
| arfken_mix_its1 | sw25 | 0.01% | 76.67% | 99.99% | 2.3 | 0.00% |
| arfken_mix_its1 | tq10 | 0.01% | 46.67% | 98.92% | 1.7 | 0.00% |
| arfken_mix_its1 | tq3 | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | tq4 | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | tq5 | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | tq6 | 0.01% | 76.67% | 98.53% | 2.3 | 0.00% |
| arfken_mix_its1 | tq8 | 0.01% | 76.67% | 99.99% | 2.3 | 0.00% |
| arfken_mix_its2 | eer | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | std | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | sw10 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | sw15 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | sw20 | 0.00% | 95.56% | 100.00% | 2.8 | 0.00% |
| arfken_mix_its2 | sw25 | 0.00% | 95.56% | 100.00% | 2.8 | 0.00% |
| arfken_mix_its2 | tq10 | 0.00% | 97.78% | 100.00% | 2.9 | 0.00% |
| arfken_mix_its2 | tq3 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | tq4 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | tq5 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | tq6 | 0.00% | 100.00% | 100.00% | 3.0 | 0.00% |
| arfken_mix_its2 | tq8 | 0.00% | 97.78% | 100.00% | 2.9 | 0.00% |
| atcc_its1 | eer | 0.00% | 88.89% | 37.13% | 5.3 | 49.74% |
| atcc_its1 | std | 0.00% | 83.33% | 37.03% | 5.0 | 50.30% |
| atcc_its1 | sw10 | 0.00% | 94.44% | 36.75% | 5.7 | 49.90% |
| atcc_its1 | sw15 | 0.00% | 100.00% | 37.76% | 6.0 | 49.40% |
| atcc_its1 | sw20 | 0.00% | 94.44% | 37.33% | 5.7 | 49.11% |
| atcc_its1 | sw25 | 0.00% | 83.33% | 36.33% | 5.0 | 49.31% |
| atcc_its1 | tq10 | 0.00% | 83.33% | 36.95% | 5.0 | 48.93% |
| atcc_its1 | tq3 | 0.00% | 83.33% | 37.03% | 5.0 | 50.30% |
| atcc_its1 | tq4 | 0.00% | 83.33% | 37.03% | 5.0 | 50.30% |
| atcc_its1 | tq5 | 0.00% | 83.33% | 37.03% | 5.0 | 50.30% |
| atcc_its1 | tq6 | 0.00% | 83.33% | 37.03% | 5.0 | 50.30% |
| atcc_its1 | tq8 | 0.00% | 94.44% | 37.42% | 5.7 | 49.30% |
| atcc_its2 | eer | 0.00% | 66.67% | 100.00% | 4.0 | 48.39% |
| atcc_its2 | std | 0.00% | 66.67% | 100.00% | 4.0 | 41.74% |
| atcc_its2 | sw10 | 0.00% | 72.22% | 100.00% | 4.3 | 47.17% |
| atcc_its2 | sw15 | 0.00% | 94.44% | 100.00% | 5.7 | 52.18% |
| atcc_its2 | sw20 | 0.00% | 77.78% | 100.00% | 4.7 | 50.87% |
| atcc_its2 | sw25 | 5.99% | 66.67% | 94.01% | 4.0 | 40.50% |
| atcc_its2 | tq10 | 0.00% | 55.56% | 100.00% | 3.3 | 56.00% |
| atcc_its2 | tq3 | 0.00% | 66.67% | 100.00% | 4.0 | 41.74% |
| atcc_its2 | tq4 | 0.00% | 66.67% | 100.00% | 4.0 | 41.74% |
| atcc_its2 | tq5 | 0.00% | 66.67% | 100.00% | 4.0 | 41.74% |
| atcc_its2 | tq6 | 0.00% | 66.67% | 100.00% | 4.0 | 41.74% |
| atcc_its2 | tq8 | 0.00% | 77.78% | 100.00% | 4.7 | 50.30% |

## Verdicts

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
