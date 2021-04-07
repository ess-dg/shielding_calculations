# Shielding calculations

This repository contains tools for calculating neutron absorption and scattering in different shielding materials. It shows how the shielding performance depends on shielding thickness, incident neutron energy, and material composistion, such as B-10 content in B4C and component ratio in epoxy-Gd2O3. In addition, it contains calculations for figure-of-merit to decide appropriate shielding parameters for shielding in different regions in the detector, which include side, back and interal shielding.

## Requisties
- Python3 (https://www.python.org/downloads/)
- ESS-DG framework (https://confluence.esss.lu.se/display/DG/CodingFramework)

## Installation

Install the ESS-DG framework according to the instructions on the webpage above.

Then, clone the repository according to:
```
git clone https://github.com/ess-dg/shielding_calculations.git
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

## Basic usage
Using the IdealGasBuilder (https://confluence.esss.lu.se/display/DG/IdealGasBuilder) any type of material can be created. The average atomic cross-sections from these materials can then be extracted according to:
```
ess_g4xsectdump_query -pneutron -lQGSP_BIC_HP_EMZ -m"IdealGas:formula=0.5*C54H60O9+0.5*Gd2O3{bymass}"
```
where in this example an epoxy-Gd2O3 shielding material was used. The generated cross-sections from the ESS-DG script is then used as input for the 'shielding_calculations'-scripts. 

The atomic densities are then calculated by hand, based on table values.
