# Shielding calculations

This repository contains tools for calculating neutron absorption and scattering in different shielding materials. It shows how the shielding performance depends on shielding thickness, incident neutron energy, and material composistion (such as B-10 content in B4C and component ratio in epoxy-Gd2O3). In addition, it contains calculations for figure-of-merit to decide appropriate shielding parameters for shielding in different regions in the detector (i.e. side-, back- and interal-shielding).

## Requisties
- Python3 (https://www.python.org/downloads/)
- ESS-DG framework (https://confluence.esss.lu.se/display/DG/CodingFramework)

## Installation

Install the ESS-DG framework according to the instructions in the webpage above.

Clone the repository:
```
git clone https://github.com/ess-dg/mg_analysis_notebook.git
```

## Execution

Open a terminal. If anaconda is installed, type:
```
conda deactivate
```

Then, type:
```
. dg_dgcode/bootstrap.sh
```

Finally, navigate to shielding_calculations/scripts and run the desired script, for example:
```
python plot_absorption.py
```
