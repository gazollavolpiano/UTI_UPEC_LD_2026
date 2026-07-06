# UTI_UPEC_LD_2026

## Study overview

This repository contains the analysis code supporting the manuscript examining whether **urinary tract infection (UTI)** and **gut carriage of uropathogenic *Escherichia coli* (UPEC)** are prospectively associated with **incident liver disease (ICD-10: K70–K77)** across **three population-based cohorts** on two continents.

- **Three cohorts** analysed through complementary registry and metagenomic pipelines: the **UK Biobank** (*n* = 479,546; clinically diagnosed UTI as the exposure), **FINRISK 2002** (*n* = 7,226; discovery cohort for gut UPEC carriage) and **SIMPLER** (*n* = 6,147; independent replication cohort)
- **UTI as a clinical exposure** (ICD-10: N39.0) associated with incident liver disease in the UK Biobank via **Cox proportional-hazards** models with sequential covariate adjustment and five prespecified sensitivity analyses, and a **four-state, five-transition multi-state** survival model.
- **Gut UPEC carriage inferred from shotgun metagenomes** and resolved through several convergent lines of evidence: a genome-wide **functional gene-family scan** (P-fimbrial *pap* genes as the discovery signal), an expanded **virulence-factor (VF) repertoire** profiled against VFDB, **phylogroup** composition (B2, D, F), and **strain-level typing** 
- **Microbial source localisation** reference-genome VF profiling across **394,913 GTDB R214 genomes**, and a **1,044-genome *E. coli* reference phylogeny** spanning the **252 sequence types** detected across both metagenomic cohorts

## Analysis pipeline summary

```
────────────────────────────  UK BIOBANK (UTI → liver disease)  ────────────────────────────
Linked ICD-10 records (n = 486,748 after consent filtering)
        │
        ├──► Prevalent/incident case definition; covariate assembly          [01_UKB/01]
        │         Exposure: prior UTI (N39.0) │ Outcome: incident LD (K70–K77)
        │         Exclude prevalent LD + viral hepatitis (B15–B19) → n = 479,546
        │
        ├──► Cox proportional-hazards models (Models 0–4, sequential)         [01_UKB/02]
        │         + 5 prespecified sensitivity analyses; Schoenfeld PH checks
        │         └──► Forest plots [forrest_plot_UKB.R, plot_all_models.R]
        │
        └──► Multi-state model — 4 states / 5 transitions (mstate)            [01_UKB/03]
                  Healthy → {diabetes, liver disease, death}; n = 450,586


──────────────────────  FINRISK 2002 (DISCOVERY)  │  SIMPLER (REPLICATION)  ──────────────────────
Stool shotgun metagenomes
   FINRISK 2002: Atropos → Bowtie2 (GRCh38) → Kraken2 (HPRC) → Kraken2 (GTDB R214) → Bracken
   SIMPLER:      CHAMP 4.0 profiler (Bowtie2 host removal → BWA-MEM → HMR05 gene catalogue)
        │
        ├──► Functional scan for incident liver disease
        │         FINRISK: HUMAnN3 gene families → Cox (detection + abundance)   [02_FR02/01]
        │                  → correlation of hits with gut E. coli abundance
        │         SIMPLER: targeted pap-operon KO test (papA–papK)               [03_SIMPLER/01]
        │
        ├──► Virulence-factor (VFDB) profiling
        │         DIAMOND blastx → RPKM → per-VF richness / burden
        │         ├──► predictor construction                      [02_FR02/02_01, 03_SIMPLER/02_01]
        │         ├──► HH/LL richness–burden grouping (k-means)
        │         ├──► Cox models: E. coli, VF richness, VF burden, HH vs LL,    [02_FR02/02_02,
        │         │        joint E. coli-adjusted models                          03_SIMPLER/02_02]
        │         ├──► single-VF Cox models → 22-VF signal set     [02_FR02/04, 03_SIMPLER/04]
        │         └──► HH-vs-LL heatmaps + ZicoSeq differential abundance
        │                                                          [02_FR02/03, 03_SIMPLER/03]
        │
        └──► E. coli strain-level typing (StrainScan)
                  FINRISK: low-depth params (--strain_prob 1 --low_dep 2)        [02_FR02/06]
                  SIMPLER: default params                                        [03_SIMPLER/05]
                  └──► strain clusters → sequence types (ST) → cohort prevalence


─────────────────────  E. coli REFERENCE-GENOME ANALYSES (Baker HPC)  ─────────────────────
GTDB R214 genomes
        │
        ├──► DIAMOND VF profiling of 394,913 bacterial genomes vs VFDB setA     [04_BAKER_HPC/01_01]
        │         ├──► VF prevalence across E. coli reference genomes           [04_BAKER_HPC/01]
        │         └──► prevalence of the 22 signal VFs across reference genomes [04_BAKER_HPC/02]
        │
        ├──► E. coli typing: Achtman MLST (mlst) + Clermont phylotyping (ezclermont)
        │                                                                       [04_BAKER_HPC/00]
        └──► Representative phylogeny of the 252 STs found across both cohorts   [04_BAKER_HPC/03]
                  1,044 CheckM-selected genomes + E. fergusonii outgroup
                  → SKA2 alignment (k = 31) → RapidNJ neighbour-joining tree
                  → Microreact (ST, CC, phylogroup, VF carriage annotation)
```

## Repository structure

```
01_UKB/                             UK Biobank — clinical UTI → incident liver disease
  01_extract_data_UKB.md              ICD-10 extraction, prevalent/incident definitions, covariates
  02_cox_models_UTI_to_LD.md          Cox PH models (Models 0–4) + sensitivity analyses
  03_multi_state_models_UTI_to_LD.md  Multi-state model (mstate)
  forrest_plot_UKB.R                  Forest plot of the fully adjusted model
  plot_all_models.R                   Combined model plots

02_FR02/                            FINRISK 2002 — discovery cohort (gut UPEC carriage)
  01_gene_fam_analysis/               HUMAnN3 gene-family scan + gut E. coli correlation
  02_VF_aggregated_analysis/          VF predictor construction + aggregate Cox models
  03_heatmap_HH_LL/                   HH/LL heatmap layers + ZicoSeq differential abundance
  04_single_VFs_analysis/             Per-VF Cox models → the 22-VF signal set
  05_VF_in_GTDB/                      Taxa list present in FINRISK metagenomes
  06_strain_analysis/                 StrainScan (low-depth) → ST typing & prevalence

03_SIMPLER/                         SIMPLER — independent replication cohort
  01_pap_KO_analysis/                 Targeted pap-operon KO replication (KEGG profiles)
  02_VF_aggregated_analysis/          VF predictor construction + aggregate Cox models
  03_heatmap_HH_LL/                   HH/LL heatmap layers + ZicoSeq differential abundance
  04_single_VFs_analysis/             Targeted replication of the 22 FINRISK VFs
  05_strain_analysis/                 StrainScan (default) → ST typing & prevalence

04_BAKER_HPC/                       E. coli reference-genome analyses (Baker HPC)
  00_Escherichia_coli_typing.md       Achtman MLST + Clermont phylotyping of GTDB E. coli
  01_01_VF_GTDB.md                    DIAMOND VF profiling of 394,913 GTDB genomes
  01_VF_Ecoli_preval.md               VF prevalence across E. coli reference genomes
  02_22_VF_GTDB_preval_plots.md       Prevalence of the 22 signal VFs across reference genomes
  03_Escherichia_coli_typing_ST_tree.md  SKA2 alignment + RapidNJ phylogeny of the 252 STs

LICENSE.md                          MIT license
README.md                           This file
```

## Software and versions

| Tool | Version | Purpose |
|---|---|---|
| NCBI datasets tool | — | Retrieval of GTDB R214 reference genomes |
| Atropos | — | Read quality-trimming and adapter filtering (FINRISK 2002) |
| Bowtie2 | 2.4.2 | Host-read depletion (mapping to human reference GRCh38) |
| Kraken2 | 2.1.3 | Taxonomic classification (HPRC human depletion + custom GTDB R214 index) |
| KrakenTools | — | Extraction of non-human reads |
| Bracken | 2.9 | Species-level abundance re-estimation (read length 150, min 25 reads) |
| HUMAnN | 3 | Functional profiling — UniRef90 gene families regrouped to Pfam domains (FINRISK 2002) |
| CHAMP profiler | 4.0 | SIMPLER read processing and functional (KEGG KO) profiling (Clinical Microbiomics) |
| BWA-MEM | 0.7.17 | Read mapping within the CHAMP pipeline (SIMPLER) |
| DIAMOND | 2.1.6.160 | Virulence-factor profiling against VFDB (metagenomes and reference genomes) |
| StrainScan | 1.0.14 | *E. coli* strain-cluster / sequence-type profiling from metagenomes |
| mlst | 2.23.0 | Achtman MLST typing of *E. coli* reference genomes |
| ezclermont | 0.7.0 | *In silico* Clermont phylotyping |
| SKA2 | 0.3.7 | Split k-mer alignment of *E. coli* genomes (*k* = 31) |
| RapidNJ | 2.3.2 | Neighbour-joining phylogenetic inference |
| CheckM | — | Genome quality scoring for representative-genome selection |
| R `survival` | 3.8.6 / 3.5.5 | Cox proportional-hazards survival models |
| R `coxme` | 2.2.18.1 | Mixed-effects Cox models (random sequencing-plate intercept, SIMPLER) |
| R `mstate` | 0.3.3 | Multi-state survival modelling (UK Biobank) |
| R `GUniFrac` (ZicoSeq) | 1.8 | Differential-abundance testing between VF HH and LL groups |
| Microreact | — | Interactive visualisation of the *E. coli* reference phylogeny |

## Reference databases

| Database | Release / date | Use |
|---|---|---|
| Human Pangenome Reference Consortium (HPRC) | k2_HPRC_20230810 | Residual human-read removal (Kraken2) |
| GTDB | R214 | Custom Kraken2 taxonomic index; *E. coli* reference-genome panel |
| VFDB | setA (downloaded 11 April 2023) | Virulence-factor reference (metagenome + genome profiling) |
| PubMLST | ecoli_achtman_4 | *E. coli* MLST scheme and clonal-complex assignment |

## Interactive resource

*E. coli* reference phylogeny with sequence-type, clonal-complex, Clermont phylogroup and 22-VF carriage annotation:
<https://microreact.org/project/e-coli-reference-phylogeny-vfs>

## Citation and license

All analysis code in this repository is released under the **MIT license** (see `LICENSE.md`) and archived at Zenodo (DOI to be added).

If you use this code, please cite: (manuscript DOI to be added).
