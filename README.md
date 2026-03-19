# Deep Brain Stimulation and Neurotransmitter Changes in Treatment-Resistant Depression

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19117530.svg)](https://doi.org/10.5281/zenodo.19117530)

This repository contains data and code supporting:

**Shevlin, Fu et al. (2026)** "Dopamine and serotonin transients predict depressive symptom relief following deep brain stimulation of human subcallosal cingulate cortex"

## Overview

The dataset includes behavioral data from 10 participants who completed both tasks described in the paper: reinforcement learning (RL) and ultimatum game (UG). The dataset also includes trial-level estimates of neurotransmitter concentrations. This repository contains implementations of the statistical analyses described in the paper as well as code for generating all figures.

## Environment

This project was developed using:
- **R version:** 4.4.3 (2025-02-28 ucrt)
- **Platform:** x86_64-w64-mingw32/x64 (Windows 11 x64, build 26100)

## Setup

First, clone this repository and navigate to the root directory:

```bash
git clone https://github.com/blairshevlin/Shevlin-Fu_Volt-TRD.git
cd Shevlin-Fu_Volt-TRD
```

All required packages are managed through `setup.R`. Run this script before any analysis:

```r
source("setup.R")
```

Installation should take less than one minute.

## Data Structure

All data are in the `data/` folder. All data have been deidentified.

- `data/behavior/` - Behavioral data for each task (RL and UG)
- `data/clinical/` - Clinical assessments for each participant
- `data/nt/` - Neurotransmitter estimates at 10Hz sampling rate
  - **Note:** Until publication, only processed data (`data/nt/processed/`) will be available
- `data/figures/` - Figure source data CSVs (one per figure, see below)

## Figure Source Data

Per Nature journal transparency requirements, the data underlying each figure's graphical representations are provided as CSV files in `data/figures/`. These are the files submitted as Source Data with the manuscript. The CSVs are not tracked by git (see `.gitignore`) — run `src/R/export_sourcedata.R` to generate them locally.

Each CSV contains:
- `panel` — identifies the sub-panel within the figure (e.g., `"A"`, `"B"`, `"RL_choice"`)
- `level` — `"trial"` (individual trial data) or `"subject"` (per-participant summaries)
- Additional columns specific to each figure (neurotransmitter estimates, behavioral measures, clinical scores, etc.)

| Figure | Source data file |
|---|---|
| Figure 1 Panel C | `data/figures/figure1_panel_source_data.csv` |
| Figure 2 | `data/figures/figure2_source_data.csv` |
| Figure 3 | `data/figures/figure3_source_data.csv` |
| Figure 4 | `data/figures/figure4_source_data.csv` |
| Extended Data Figure 1 | `data/figures/figureExtended1_source_data.csv` |
| Extended Data Figure 2 | `data/figures/figureExtended2_source_data.csv` |
| Supplementary Figure 1 | `data/figures/figureSupplement1_source_data.csv` |
| Supplementary Figure 2 | `data/figures/figureSupplement2_source_data.csv` |
| Supplementary Figure 3 | `data/figures/figureSupplement3_source_data.csv` |
| Supplementary Figure 4 | `data/figures/figureSupplement4_source_data.csv` |
| Supplementary Figure 5 | `data/figures/figureSupplement5_source_data.csv` |
| Supplementary Figure 6 | `data/figures/figureSupplement6_source_data.csv` |

## Scripts

All code for generating figures and conducting statistical analyses are in the `src/` folder:

- `process_nt.r` - Converts raw estimates into processed data used for all analyses
  - **Note:** Requires access to `data/nt/raw/`. This script will not run until raw data is released upon publication
- `statistical_tests.R` - Runs all statistical analyses reported in the main text
- `export_sourcedata.R` - Exports figure source data CSVs to `data/figures/` (run once before generating figures)
- `generate_figureXX.R` - Generates the corresponding figure from its source data CSV (XX = figure number from main text). Output figures are saved to `results/from_source_data/`.

Run time for each script is less than 5 minutes.

## Usage

1. Run `source("setup.R")` to install and load required packages
2. Export figure source data: `source("src/R/export_sourcedata.R")` — writes all CSVs to `data/figures/`
3. Execute figure generation scripts: `source("src/R/generate_figure2.R")`, etc. — output saved to `results/from_source_data/`
4. Run statistical analyses: `source("src/R/statistical_tests.R")`

## Citation

Individuals may use or adapt the code provided they follow the terms of the license. When using the code, please cite:

### BibTeX:
```bibtex
@software{shevlin2026analysis,
  author = {Shevlin, Blair and Fu, Qi Xiu},
  title = {Analysis code for dopamine and serotonin transients predict depressive symptom relief following deep brain stimulation of human subcallosal cingulate cortex},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.19117530},
  url = {https://doi.org/10.5281/zenodo.19117530}
}
```

### APA:
```
Shevlin, B., & Fu, Q. X. (2026). Analysis code for dopamine and serotonin transients predict depressive symptom relief following deep brain stimulation of human subcallosal cingulate cortex. Zenodo. https://doi.org/10.5281/zenodo.19117530
```

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.