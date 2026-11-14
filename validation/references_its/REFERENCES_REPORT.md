# Reference ITS sequences from NCBI

In-silico PCR with the primers of Arfken et al. 2023 (Table 1); products without the primers. 'strain': an assembly of the very strain; 'species': another strain of the species; 'given': chosen with --use.

| Member | Mock | Species | Assembly | Strain in the assembly | ITS1 reference is | ITS2 reference is | ITS1 products (nt x copies) | ITS2 products (nt x copies) | Edits to the isolate (ITS1 / ITS2) | Source |
|---|---|---|---|---|---|---|---|---|---|---|
| Af | ATCC | Aspergillus fumigatus | GCF_000002655.1 | Af293 | strain (genome) | strain (genome) | 260x1 | 313x1 |  | genome GCF_000002655.1 |
| Cn | ATCC | Cryptococcus neoformans var. grubii | GCF_000149245.1 | H99 | strain (genome) | strain (genome) | 197x1 | 334x1 |  | genome GCF_000149245.1 |
| Ti | ATCC | Trichophyton interdigitale | GCA_037576225.1 | IGIB-BK-CI04 | species (genome) | species (genome) | 327x2, 324x1, 325x1, 321x1 | 322x1, 317x1, 320x1, 319x1, 322x1 |  | genome GCA_037576225.1 |
| Pc | ATCC | Penicillium chrysogenum | GCF_028827035.1 | IBT 35668 | species (genome) | species (genome) | 251x1 | 311x1 |  | genome GCF_028827035.1 |
| Fk | ATCC | Fusarium keratoplasticum | GCF_025433545.1 | Fu6.1 | species (genome) | species (genome) | 227x7, 200x1, 229x1, 241x1, 231x1, 229x1, 232x1, 230x1, 230x1, 229x1, 232x1, 230x1 | 318x4, 318x1, 318x1, 320x1, 319x1, 309x1, 321x1, 320x1, 321x1, 322x1, 319x1, 320x1, 322x1, 321x1, 317x1 |  | genome GCF_025433545.1 |
| Ca | ATCC | Candida albicans | GCA_002276455.1 | ATCC 10231 | strain (record) | strain (record) | 214x1 | 298x3 |  | genome GCA_002276455.1; ITS1 from NCBI ITS records of 10231; ITS2 from NCBI ITS records of 10231 |
| Cg | ATCC | Nakaseomyces glabratus | GCF_010111755.1 | ATCC 2001 | strain (genome) | strain (genome) | 478x6, 477x2 | 377x8 |  | genome GCF_010111755.1 |
| Mg | ATCC | Malassezia globosa | GCF_000181695.2 | CBS 7966 | strain (genome) | strain (genome) | 328x1 | 436x1 |  | genome GCF_000181695.2 |
| Sc | ATCC | Saccharomyces cerevisiae | GCF_000146045.2 | S288C | strain (genome) | strain (genome) | 437x2 | 380x2 |  | genome GCF_000146045.2 |
| Cd | ATCC | Cutaneotrichosporon dermatis | GCA_020011015.1 | NICC30027 | species (genome) | species (genome) | 193x3 | 311x3 |  | genome GCA_020011015.1 |
| C | piglet | Candida albicans | GCF_000182965.3 | SC5314 | species (genome) | species (genome) | 214x1 | 298x1 | 1 (isolate 215 nt) / 0 (isolate 298 nt) | genome GCF_000182965.3 |
| K | piglet | Kazachstania slooffiae | GCA_017347545.2 | ABBL | species (genome) | species (genome) | 395x6, 394x1 | 390x6, 390x1 | 0 (isolate 395 nt) / 0 (isolate 390 nt) | genome GCA_017347545.2 |
| M | piglet | Lichtheimia corymbifera | GCA_000697475.1 | B2541 | none | species (record) | none | 393x1 | - / 56 (isolate 399 nt) | genome GCA_000697475.1; ITS2 from NCBI ITS records of the species |
| P | piglet | Pichia fermentans | GCA_035610515.1 | CDLB.YE03 | species (genome) | species (genome) | 154x17, 154x11 | 266x10, 267x8, 267x6, 266x1, 267x1 | 0 (isolate 154 nt) / 3 (isolate 266 nt) | genome GCA_035610515.1 |
| T | piglet | Trichosporon asahii | GCF_000293215.1 | CBS 2479 | species (genome) | species (genome) | 199x1 | 318x1 | 0 (isolate 199 nt) / 1 (isolate 318 nt) | genome GCF_000293215.1 |

Members without any product: none.

ATCC members with a reference of the very strain (from its genome or its own ITS records): ITS1 6 of 10 (Af, Cn, Ca, Cg, Mg, Sc); ITS2 6 of 10 (Af, Cn, Ca, Cg, Mg, Sc).

Read before the plan uses these: a member with no product, a product far from the length the paper or run 25 shows, or a piglet member many edits from its isolate needs another assembly (--use) or stays out of the truth. Next: send this file.
