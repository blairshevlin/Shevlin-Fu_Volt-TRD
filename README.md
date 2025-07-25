# Deep Brain Stimulation and Neurotransmitter Changes in Treatment-Resistant Depression

This repository contains data and code supporting:

**Shevlin, Fu et al. (2025)** "Deep brain stimulation to the subcallosal cingulate induces context-dependent changes in dopamine and serotonin in humans with treatment-resistant depression"

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

## Scripts

All code for generating figures and conducting statistical analyses are in the `src/` folder:

- `process_nt.r` - Converts raw estimates into processed data used for all analyses
  - **Note:** Requires access to `data/nt/raw/`. This script will not run until raw data is released upon publication
- `statistical_tests.R` - Runs all statistical analyses reported in the main text
- `generate_figureXX.R` - Generates the corresponding figure (XX = figure number from main text).

Run time for each script is less than 5 minutes.

## Usage

1. Run `source("setup.R")` to install and load required packages
2. Execute figure generation scripts: `source("generate_figure2.R")`, etc.
3. Run statistical analyses: `source("statistical_tests.R")`

## Citation

Individuals may use or adapt the code provided they follow the terms of the license. When using the code, please cite:

### BibTeX:
```bibtex
@misc{shevlin2025analysis,
  author = {Shevlin, Blair and Fu, Qi Xiu},
  title = {Analysis code for deep brain stimulation to the subcallosal cingulate induces context-dependent changes in dopamine and serotonin in humans with treatment-resistant depression},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/blairshevlin/Shevlin-Fu_Volt-TRD}
}
```

### APA:
```
Shevlin, B., & Fu, Q. X. (2025). Analysis code for deep brain stimulation to the subcallosal cingulate induces context-dependent changes in dopamine and serotonin in humans with treatment-resistant depression. GitHub. https://github.com/blairshevlin/Shevlin-Fu_Volt-TRD
```

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.