# Neural Dynamics of Abstractness and Abstraction in Word Processing

This repository contains the analysis code and data accompanying the paper:

> **Unraveling the interplay between sensory information and semantic architecture in brain concept representation: the role of abstraction**
> **Authors: Chiara Battaglini, Davide Bottari, Giacomo Handjaras, Martina Berto, Ella Striem-Amit,**
> **Pietro Pietrini, Alessandro Lenci, Giovanna Marotta, Emiliano Ricciardi**

---

## Overview

Understanding how abstract and concrete words are represented in the brain requires disentangling the contribution of sensorimotor information from the organization of semantic knowledge. This project introduces a distinction between **abstractness** (sensorimotor grounding) and **abstraction** (semantic architecture), which are typically confounded in the literature. By including superordinate words as a critical comparison category, we dissociate these dimensions and identify their distinct neural correlates using EEG and Representational Similarity Analysis (RSA).

Key findings:
- Superordinate words pattern with concrete words on sensory properties but with abstract words on measures of semantic architecture.
- Concrete words elicit stronger negativities in the N400 and N700 time windows compared to both abstract and superordinate words.
- RSA reveals that neural response patterns are best explained by an abstraction model, outperforming models based on sensorimotor or distributional information.

---

## Repository Structure

```
CODE/
├── functions/
│   ├── compute_corr.m                    # Correlation computation utilities
│   ├── get_sig_windows.m                 # Identifies significant time windows
│   └── plot_cluster_markers.m            # Plots cluster-based significance markers
│   └──fdr_bh.m                           # False discovery rate correction
├── ERP_grandaverage_plot.m               # Grand average ERP visualization
├── cluster_permutation_and_plot.m        # Cluster-based permutation test + plotting
├── cluster_permutation_statistics_3...m  # Omnibus cluster-based permutation test
└── timeresolved_RSA.m                    # Time-resolved Representational Similarity Analysis

Data/
├── Behavioral/
│   ├── EEG_subsample.csv                 # Behavioral data for the EEG subsample
│   └── Whole_database.csv                # Full behavioral norming database
├── EEG/                                  # averaged EEG data
└── RSA/                                  # RSA model matrices and results
```

---

## Methods

### Behavioral Norming
Feature production and norming tasks (abstractness, concreteness, familiarity, generalizability) were used to quantify sensory-related properties and semantic architecture across word categories (concrete, abstract, superordinate). From feature production feature relevance, and positive pointwise mutual information (PPMI) were extracted.

### EEG Acquisition & Preprocessing
Neural data were acquired using electroencephalography (EEG) with an EGI's HydroCel Geodesic Sensor Net with 65 EEG channels and a Net Amps 400 amplifier (Electrical Geodesics, Inc., EGI, USA) while participants heard words. Standard semiautomatic preprocessing pipelines were applied (details in the paper).

### Representational Similarity Analysis (RSA)
Time-resolved RSA was used to assess the correspondence between computational models and neural response patterns over time. Models capture:
- **Abstractness** – concreteness and abstactness of words
- **Abstraction** – semantic architecture (generalizability, relevance and PPMI)
- **Distributional similarity** – word co-occurrence statistics in the ItWac128 corpus using Word2Vec 
- **Distributional similarity** – sensorimotor norms from The Lancaster Sensorimotor Norms

---

## Requirements

- **MATLAB** (developed and tested on R2021b or later)
- **EEGLAB** toolbox (for EEG data handling)
- **FieldTrip** toolbox (for cluster-based permutation statistics)

No additional toolboxes are required beyond those listed above.

---

## How to Use

### 1. Clone the repository

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

### 2. Set up MATLAB paths

In MATLAB, add the repository and required toolboxes to your path:

```matlab
addpath(genpath('CODE/'))
addpath(genpath('/path/to/eeglab'))
addpath(genpath('/path/to/fieldtrip'))
```

### 3. Run ERP analyses

To reproduce the grand average ERP plots:

```matlab
ERP_grandaverage_plot
```

To run the cluster-based permutation statistics and generate ERP comparison plots:

```matlab
cluster_permutation_and_plot
```

### 4. Run time-resolved RSA

To compute model–brain correlations over time and obtain RSA results:

```matlab
timeresolved_RSA
```

### 5. Explore behavioral data

The `Data/Behavioral/` folder contains:
- `EEG_subsample.csv` — median of behavioral ratings and feature measures for stimuli used in the EEG study
- `Whole_database.csv` — the full norming database including all stimuli

These can be loaded in MATLAB or any standard data analysis tool (R, Python, Excel).


---

## Citation

If you use this code or data in your work, please cite:

```
[Author(s)] (Year). Title of the paper. Journal Name, Volume(Issue), Pages.
DOI: https://doi.org/xxxx
```

> *(Citation will be updated upon publication.)*

---

## License

This code is released for academic use. Please contact the authors for other use cases.

---

## Contact

For questions about the code or data, please open an issue or contact the corresponding author at chiara.battaglini@iusspavia.it.
