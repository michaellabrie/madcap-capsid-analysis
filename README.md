# Sequence-Based Discovery and Prioritization of Novel AAV Capsids

## Project Overview

This solo project analyzes capsid sequence-identity data from the 2026 MADCAP study by Lyashenko et al.

The dataset compares 139 MADCAP-discovered adeno-associated virus (AAV) capsids with 162 previously known AAV isolates. The goal is to examine sequence diversity and identify unusual capsids that may be useful candidates for future structural and experimental characterization.

This project was completed for CAI5735: Applied Data Science in Health at the University of Florida.

## Research Question

How are nearest-known sequence identity and overall sequence similarity related among MADCAP-discovered capsids, and does this relationship differ between divergent and similar capsids?

## Dataset

Primary dataset:

MADCAP Supplementary Table S5 - Identity matrix of all capsids used in the study.

Source:

Lyashenko et al. (2026), Journal of Virology.

"Metagenomic AAV Discovery and Characterization Pipeline (MADCAP): isolation of novel nAb-naive AAV capsids from metagenomic data."

The dataset contains:

- 301 total capsids
- 139 MADCAP-discovered capsids
- 162 previously known isolates
- 302 columns in the source CSV

## Analysis

The Marimo notebook includes:

- filtering MADCAP capsids from known isolates
- sorting capsids by sequence divergence
- grouping capsids by divergence status
- computed sequence-diversity columns
- histogram of nearest-known sequence identity
- boxplot comparing divergence groups
- scatterplot of nearest-known identity vs. mean-known identity
- group comparison of divergent and similar capsids
- ranking of the most divergent candidates

## Key Results

- 63 MADCAP capsids were classified as divergent using the study's <95% identity threshold.
- 76 capsids had >=95% identity to their nearest known isolate.
- Nearest-known sequence identity ranged from approximately 58.4% to 99.9%.
- Nearest-known identity and average identity across known isolates were strongly related (Pearson r approximately 0.79).
- MC68 was the most sequence-divergent candidate in the analysis.

Sequence divergence is used here as a prioritization metric, not as a prediction of biological performance.

## Repository Structure

- `MADCAP_Solo_Marimo_Notebook.py` - final solo Marimo analysis
- `madcap_notebook_original_backup.py` - backup of the earlier notebook
- `data/raw/` - original MADCAP supplementary dataset
- `figures/` - analysis figures
- `requirements.txt` - Python dependencies

## Running the Notebook

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
