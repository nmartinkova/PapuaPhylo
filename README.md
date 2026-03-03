# Elevational niche evolution in New Guinean murine rodents

This repository contains R scripts and associated files used in the analyses for:

**Directional dispersal from an ancient montane core drives diversification across New Guinea’s mountain systems**

The scripts reproduce the phylogenetic, ancestral state, and biogeographical analyses presented in the manuscript.

---

## Repository structure

* `1map.r` -- Generates sampling maps and visualisations of elevational transects in New Guinea.
* `2clockCheck.r` -- Performs molecular clock diagnostics and model comparison prior to time calibration.
* `3ancestralStates.r` -- Implements ancestral elevation reconstructions and discrete regional state analyses, including stochastic character mapping.
* `4stratifiedSubsampling.r` -- Conducts stratified subsampling analyses to assess robustness to unequal lineage representation between clades.
* `5simmapStrata.tex` -- LaTeX code used to format and visualise stochastic mapping summaries in the manuscript.

---

## Data

Raw genomic data (ddRAD-derived SNP datasets) are deposited in [repository name / accession number – to be inserted].

Processed phylogenetic trees and intermediate objects required to run the scripts should be placed in the working directory as specified within each script.

---

## Software and packages

Analyses were conducted in **R (4.5.2)**.
Key packages include:

* `ape`
* `phytools`
* `OUwie`

Additional package dependencies are listed within each script.

---

