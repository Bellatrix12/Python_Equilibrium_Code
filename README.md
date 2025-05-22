# Python_Equilibrium_Code
### Overview

This set of python modules and databases calculate the equilibrium state for a given planetary system, and returns the Gibbs free energy difference
between the initial (observed) state and the equilibrium states. Please see J. Krissansen-Totton, S. Olson and D.C. Catling (2018) "Disequilibrium biosignatures over Earth history and implications for detecting exoplanet life", Science Advances, 4, eaao5747 and J. Krissansen-Totton, D. S. Bergsman and D.C. Catling (2016) "On Detecting Biospheres from Thermodynamic Disequilibrium in Planetary Atmospheres", Astrobiology 16, 39-67 for a detailed description of the science behind these code modules.

These python modules were translated and adapted by Amber V. Young using the original Matlab scripts developed by Krissansen-Totton et al. We translated the original Matlab scripts so that they can be run on open source software. **Note that at this time, the python version can
only run gas phase systems.**

The examples provided are for Mars and Modern Earth, which provide results that are consistent with results published in Krissansen-Totton et al., 2016.

### Requirements
Python 3.7 or later (earlier versions of python may work but are untested)

### How to Run the Code
      (1) Make sure you are in the top level directory of the model and ensure python is installed and operational.
      (2) Run the model in one of the two following ways:
                - In a command terminal execute the command: python Main_Script_iterate.py
                - In a python IDE of your choice run the Main_Script_iterate.py module using the GUI


Example printout for modern Earth gas phase system:

```
Optimization terminated successfully    (Exit mode 0)
            Current function value: -7987711.572338882
            Iterations: 4
            Function evaluations: 68
            Gradient evaluations: 4
Mass Balance Roughness (should be zero):  0.0
dG value:  -6316.862378405407
G1 value:  -7981394.7099604765
G2 value:  -7987711.572338882
deltaG value this iteration (J/mol):  -1.5134085279469294
```


The available Gibbs free energy is printed as the deltaG value and is in units of joules per mole of atmosphere. There will also be a bar graph output
comparing the initial and equilibrium states for all chemical species included in the calculation.

To switch between atmospheric input templates make sure to read through the `load_vectors.py` module and adapt it accordingly with the appropriate file name for the inputs, and update the number of species for the variable `n_spec`. 


**For detailed descriptions of each module and function see the associated `ReadMe_AVY.txt` text file.**
