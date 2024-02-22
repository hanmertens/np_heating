# Modeling of the heating of layered plasmonic nanoparticles

## Authors

Johannes C. J. Mertens and Monique A. van der Veen; Department of Chemical Engineering, Delft University of Technology.

## Description

Code for modeling the heat transfer and calculating the temperature evolution of a layered plasmonic nanoparticle and its surroundings after pulsed or continuous-wave illumination, based on the work by Metwally *et al.*[1]
Optical constants of gold were sourced from [2], material properties from [3] and interfacial conductances from [4-8].

## Usage

The following needs to be installed in order to use this model:
- A somewhat recent version of Python 3; the model was tested using Python 3.10.
- Several Python dependencies listed in `requirements.txt`, which can be installed using `pip install -r requirements.txt`.

With these installed, the various scripts named `script_*.py` can be run.
These scripts have various parameters that can be changed near the top of their file, for example to change the core radius, shell thickness, or medium material.

## List of files

Core files:
- `np_heating.py`: Module performing heat transfer modeling.
- `nAu.txt`: Optical constants of gold,[2] used by `np_heating.py`. Space-separated values, with the first column the energy in electronvolt, with the second and third columns the dimensionless optical constants *n* and *k*, respectively.

Scripts using `np_heating.py` to model various situations:
- `script_bare_vs_shell.py`: Calculates and plots temporal and spatial temperature profiles of gold nanoparticles with and without a silica shell after femtosecond-pulsed illumination.
- `script_shell_thickness.py`: Calculates and plots temporal temperature profiles of the environment of gold nanoparticles with different silica shell thicknesses after femtosecond-pulsed illumination.
- `script_interfacial_conductance.py`: Calculates and plots temporal temperature profiles of gold nanoparticles with or without silica shell with different interfacial conductances after femtosecond-pulsed illumination.

Documentation:
- `README.md`: This file.
- `requirements.txt`: List of dependencies of `np_heating.py`, which can be used by a tool such as `pip`.

## References

1. K. Metwally, S. Mensah, G. Baffou, *J. Phys. Chem. C* **2015**, *119*, 28586.
2. P. B. Johnson, R. W. Christy, *Phys. Rev. B* **1972**, *6*, 4370.
3. W. M. Haynes (Editor), *CRC Handbook of Chemistry and Physics*, CRC Press, Boca Raton, Fla., 96th edition **2015**.
4. S. M. Hatam-Lee, F. Jabbari, A. Rajabpour, *Nanoscale Microscale Thermophys. Eng.* **2022**, *26*, 40.
5. A. Plech, V. Kotaidis, S. Grésillon, C. Dahmen, G. von Plessen, *Phys. Rev. B* **2004**, *70*, 195423.
6. J. Park, J. Huang, W. Wang, C. J. Murphy, D. G. Cahill, *J. Phys. Chem. C* **2012**, *116*, 26335.
7. P. A. Schoen, B. Michel, A. Curioni, D. Poulikakos, *Chem. Phys. Lett.* **2009**, *476*, 271.
8. J. Park, D. G. Cahill, *J. Phys. Chem. C* **2016**, *120*, 2814.
