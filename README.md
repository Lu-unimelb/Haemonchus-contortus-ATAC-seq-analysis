# ATAC-seq Analysis of Chromatin Accessibility in *Haemonchus contortus*

This repository contains the complete analysis workflow, including scripts, intermediate files, and tutorial-style notes for the investigation of chromatin accessibility in the parasitic nematode *Haemonchus contortus*.

> ** Associated Manuscript:** *placeholder for your paper's title*
>
> **DOI:** `https://doi.org/10.xxxx/xxxxx` (link to be added upon publication)

-----

## Repository Structure

  * `data/`: Directory for genome files, annotation references, and sample sheets (not included).
  * `scripts/`: Contains executable `.R` scripts and explanatory `.md` walkthroughs.
  * `results/`: All output files generated during the analysis (e.g., peaks, filtered regions, functional analysis plots).
  * `README.md`: Project overview and instructions (this file).

-----

## Data Availability

Due to their size, raw data and large reference files are **not included** in this repository. To replicate the analysis, please download the required genome and annotation files and place them in a local `data/` directory.

> **Note:** The exact versions and sources for all reference files used in the published analysis are detailed in the manuscript's Materials and Methods section.

-----

## Analysis Workflow

This project follows a "tutorial + scripted record" style. Follow the numbered scripts and markdown files in the `scripts/` directory in order. Each step corresponds to a specific section in the manuscript's **Materials and Methods**.

1.  **01 – File Preparation** (`Methods 2.2`)

      * `scripts/01_file_prep.md`
      * This markdown guide covers the preparation of all required input files, including raw sequencing data and reference genome files, for the subsequent pipeline steps.

2.  **02 – ATAC-seq Data Processing** (`Methods 2.3`)

      * `scripts/02_ATACseq_data_analysis/02_ATACseq_data_analysis.md`
      * Processes the raw ATAC-seq data using the **nf-core/atacseq** pipeline ([Ewels et al., 2020](https://doi.org/10.1038/s41587-020-0439-x)).

3.  **03 – Peak QC and Filtering** (`Methods 2.4`)

      * `scripts/03_peaks_qc_and_filtering/03_peaks_qc_and_filtering.md`
      * Performs quality control on the peaks called by the pipeline and uses the Irreproducible Discovery Rate (IDR) framework to identify a high-confidence set of reproducible peaks for downstream analysis.

4.  **04 – Defining Promoter-Accessible Genes** (`Methods 2.5`)

      * `scripts/04_promoter_accessible_genes/04_promoter_accessible_genes.md`: Defines promoter-accessible genes based on peak locations.
      * `scripts/04_promoter_accessible_genes/04_promoter_accessible_genes.R`: Performs statistical analysis and visualization of promoter-accessible genes.

5.  **05 – Functional Analysis & Enrichment** (`Methods 2.5 & 2.6`)

      * `scripts/05_functional_analysis/05_functional_analysis.md`: Generates necessary input files for enrichment and performs functional annotation using **eggNOG-mapper** ([Cantalapiedra et al., 2021](https://doi.org/10.1093/molbev/msab293)).
      * `scripts/05_functional_analysis/accessibility_and_transcription.R`: Analyzes the relationship between gene promoter accessibility and transcription levels.
      * `scripts/05_functional_analysis/enrichment_and_plot.R`: Conducts functional enrichment analysis (e.g., GO, KEGG) and generates plots.

6.  **06 – Gene Essentiality Analysis** (`Methods 2.6`)

      * `scripts/06_gene_essentiality_analysis/06_gene_essentiality_analysis.R`
      * Investigates the relationship between chromatin accessibility at gene promoters and gene essentiality scores.

7.  **07 – Accessibility in Co-expression Modules** (`Methods 2.8`)

      * `scripts/07_accessibility_in_coexpression_modules/07_accessibility_in_coexpression_modules.R`
      * Analyzes and compares chromatin accessibility across distinct gene co-expression modules.

8.  **08 – Motif Analysis** (`Methods 2.7`)

      * `scripts/08_motif_analysis/08_motif_analysis.R`
      * Performs *de novo* motif discovery and known motif enrichment analysis within the set of reproducible ATAC-seq peaks to identify potential transcription factor binding sites.

-----

## Citation

If you use the code, analysis workflow, or data from this repository in your research, please cite our manuscript. You can also cite this repository directly using the information below.

```bibtex
@misc{Lu2025_HcATAC,
  author       = {Lu, Hanyu},
  title        = {Haemonchus-contortus-ATAC-seq-analysis: Tutorial-style notes and scripts},
  year         = {2025},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/Lu-unimelb/Haemonchus-contortus-ATAC-seq-analysis}},
  note         = {See associated manuscript for methods and software versions},
  doi          = {10.xxxx/xxxxx} % <-- Replace with your paper's DOI
}
```

-----

## Contact & License

  * **Contact:** For questions, bug reports, or discussions, please open a GitHub Issue.
  * **License:** The code in this repository is available under the [MIT License](LICENSE.md). The documentation and tutorial notes are licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).