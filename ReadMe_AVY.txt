This set of python modules and databases calculate the equilibrium
state for a given planetary system, and returns the Gibbs free energy difference
between the initial (observed) state and the equilibrium states. Please see
J. Krissansen-Totton, S. Olson and D.C. Catling (2018)
"Disequilibrium biosignatures over Earth history and implications for detecting
exoplanet life", Science Advances, 4, eaao5747 and
J. Krissansen-Totton, D. S. Bergsman and D.C. Catling (2016) "On Detecting
Biospheres from Thermodynamic Disequilibrium in Planetary Atmospheres",
Astrobiology 16, 39-67 for a detailed description of the science behind these
code modules.

These python modules were translated and adapted by Amber V. Young using
the original Matlab scripts developed by Krissansen-Totton et al. .
We translated the original Matlab scripts so that they can
be run on open source software. Note that at this time, the python version can
only run gas phase systems.

The examples provided are for Mars and Modern Earth. The Mars and modern
Earth examples are consistent with results published in
Krissansen-Totton et al., 2016.

This work was supported by the NASA Exobiology Program, and NASA's Nexus for
Exoplanet Science Virtual Planetary Laboratory.

REQUIREMENTS: Python 3.7 or later (earlier versions of python may work but
are untested)
#~*~*~*~**~*~**~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

HOW TO RUN THE CODE:
(1) Make sure you are in the top level directory of the model and ensure python
is installed and operational.
(2) Run the model in one of the two following ways:
      - In a command terminal execute the command: python Main_Script_iterate.py
      - In a python IDE of your choice run the Main_Script_iterate.py module
      using the GUI


Example printout for modern Earth gas phase system:

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


The available Gibbs free energy is printed as the deltaG value and is in units
of joules per mole of atmosphere. There will also be a bar graph output
comparing the initial and equilibrium states for all chemical species included
in the calculation.

To switch between atmospheric input templates make sure to read through the
load_vectors module and adapt it accordingly with the appropriate file name
for the inputs, and update the number of species for the variable n_spec. 
#~*~*~**~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~*~*


EXPLANATION OF CODE STRUCTURE:

Main_script_iterate:
This module loads the input file (see below) and performs the global
optimization routine by calling the Gibbs energy minimization module for
multiple, randomized initial conditions. The minimum value (largest negative
change in Gibbs energy) is selected and displayed, and the plotting script
(Plot_outputs) is called. The number of iterations is specified in this module
by the variable num.

Gibbs_energy_minimization:
This script calculates the equilibrium state given some initial
(observed) out-of-equilibrium state, and calculates the Gibbs free energy
difference between the intial and equilibrium states. The following inputs are
required (the input section is clearly outlined in the code):
  - Pressure of the system (default P=1.013 bar)
  - Temperature of the system (default T=288.15 K)
  - Arrays containing the initial state species abundances (mol/kg). These are
  obtained from the load_vectors module.
  - Random initial state vectors to be used as the arbitrary starting point for
  the optimization routine. These are automatically provided by the
  Main_script_iterate module.


Atom conservation conditions will be calculated from the true initial
state. Specifically, the "makeAs" module is used to create a dictionary
comprised of the element(s) that make up a particular species.
The name of the element is stored as the dictionary key(s) and the number of
 atoms of each element is stored as a list tied to the dictionary key. By
computing the Gibbs energy of each species, the total Gibbs free energy of the
true initial (observed) state is calculated by calling the function
totalGibbsInternal (see below).

Next, the random initial condition along with the atom constraints
from the initial state are fed into the scipy.optimize.minimize routine. This
routine uses an interior points method of optimization identical to the
one used in the Matlab version of this model. The scipy.optimize.minimize
routine minimizes the total Gibbs energy of the system, as defined by the
the function totalGibbsInternal, subject to the constraint that atoms
are conserved and that mixing ratios for all species are non-negative.
In fact, we bound the mixing ratios between 1e-30 and 1 to avoid singularities
from taking the natural log of zero. The minimization routine returns a
dictionary of output labeled as the "outputs" variable in the module. We save
from those outputs the molar abundances in equilibrium and the Gibbs free
energy of the final state.

Finally given the Gibbs free energy of the final state, the difference in
Gibbs energies between the initial and final state is calculated and
printed.

Within Gibbs_energy_minimization, there are three local functions that are
called:
  - totalGibbsInternal: This takes the molar abundances array and calculates
  the Gibbs energy of the system according to equation 10 of Krissansen-Totton
  et al. (2018). This function calls various external modules like fugCoef
  (calculates the fugacity coefficients), fugCoefAQ (Calculates the activity
  coefficients for aqueous species), and Pitzer_activity_diseq (calculates
  aqueous activity coefficients using the Pitzer equations). The function
  also calculates the analytic gradient of the Gibbs energy with respect
  to abundances, which ensures the optimization routine converges more
  reliably. Again, note that this version of the model only runs for gas phase
  systems meaning any modules involving aqueous phase species go unused.

  The total Gibbs energy of the system is returned (the fun variable).

  - onlyPhase: This function returns a form of some original array with only
  species in it of the specified phase.
  -mycon: Used for calculating charge balance (goes unused).
  -mass_constraint: Used for calculating mass balance.
  -const: Used to specify all the constraints of the system


Plot_ouputs:
This module takes the initial and final abundances as inputs and plots them
as bar charts.


gibbs:
The gibbs module takes the following inputs:
  - List of names of gaseous molecular species (or the name of a single species)
  - Temperature
It computes and returns standard Gibbs free energies of formation for gaseous
species by calling the  NEWNASA database. The NEWNASA database actually
contains coefficients for computing the Berman-Brown Gibbs free energies
but the gibbs module converts these to standard Gibbs free energies of formation.
The gibbs module calls external modules such as makeAs and convert2standard.
The former creates a dictionary of each species with the elements its composed
of and the number of atoms of each element. The latter module converts
certain elements to their standard forms for the purposes of Gibbs free
energy calculations (see below).

gibbsBB:
Identical to the gibbs module except the Berman-Brown Gibbs free energies
are retained.

gibbsAQ:
This module takes the following inputs:
- list of names of aqueous molecular species (or the name of a single species)
- Temperature
- Pressure
It computes and returns Gibbs free energies of formation for aqueous species by
calling the SUPCRT database (see loadDatabaseC). See Appendix C in
Krissansen-Totton et al. 2016 for a complete description of how Gibbs energies
are calculated for aqueous species. The module goes unused for gas phase
calculations

makeAs:
The makeAs takes a cell array of species names and creates the
coefficient matrix (in the form of a dictionary) used to define the atom
conservation constraint. The name of the element is stored as the
dictionary key(s) and the number of atoms of each element is stored as a list
tied to the dictionary key.

Parsing:
The Parsing module takes in species names as inputs and determines their
constituent elements. The output is a list of elements within that species

convert2standard:
The "convert2standard" function ensures that molecular species are converted to
their reference state for the purposes of Gibbs free energy calculations
(e.g. O present as O2, carbon is present as graphite). It takes an array of
molecular species as an input, and returns the same species in their standard form.

fugCoef:
This module returns fugacity coefficient values using the Soave equation
(see Appendix A in Krissansen-Totton et al. 2016 for a complete description).
It requires the following inputs:
- Temperature and pressure of the system
- List of names of molecular species for which fugacity coefficient is to be calculated.
- Array of molar abundances of these species
The module calls databaseB which is created by the module "loadDatabaseB",
described below. databaseB contains the critical temperatures,
critical pressures and acentric factors required to calculate
fugacity coefficients. The "fugCoef" function returns a cell array containing
log(fugacity coefficients) for each of the inputted species. Errors will be
returned if the designated species cannot be located in the database or if the
computer has trouble solving the cubic equation for the compressibility factor.

fugCoefAQ:
This module calculates aqueous activity coefficients using the
Truesdell-Jones equation (see Krissansen-Totton et al. 2016 for details).
This method has been largely superseded by calculating activity coefficients
using the more accurate Pitzer equations (see below). This module goes unused
for the gas phase calculations

Pitzer_activity_diseq:
This module calculates aqueous activity coefficients and water activity using
the Pitzer equations. The methodology for calculating water activity is
described in Krissansen-Totton et al. 2016, whereas the methodology for
calculating activity coefficients for other aqueous species is described
in the 2018 paper (equation 11). In short, this
function takes as inputs the abundances, charges, and names of aqueous species
and outputs their activity coefficients. Activity coefficients for anions and
cations are calculated separately using the database of P
itzer coefficients "pitzer_coef.txt" (see load_Pitzer). The temperature
dependence of aqueous activity coefficients is neglected. This function is
called by totalGibbsInternal when calculating the total Gibbs energy of a system.
This module goes unused for the gas phase calculations.

dielectric:
Calculates the dielectric constant of water, which is used to calculate the
Gibbs energies of aqueous species in gibbsAQ.

load_vectors:
This module opens the designated input file containing species names,
abundances, charges, and phase. The input file is read and each line is
converted to an array or cell, as needed.

loadDatabase:
The "loadDatabase" module opens the NEWNASA.txt file that contains coefficients
for calculating Gibbs free energies of formation. The coefficients for all
species are stored as a Matlab array called "database" which is called by the
"gibbs" function described above.

loadDatabaseB:
The "loadDatabaseB" function opens the fugacityCoefficientVariables.txt file
that contains thermodynamic variables for calculating gas fugacity coefficients.
The variables for all species are stored as a Matlab array called "databaseB"
which is called by the "fugCoef" function described above.

loadDatabaseC:
This module opens the sprons96_edited2.dat file that contains coefficients for
calculating aqueous Gibbs free energies of formation. The variables for all
species are stored as the Matlab array called �databaseC� which is called by
the �gibbsAQ� function described above.

loadDatabaseD:
This module opens the AQ_coefficients.txt file that contains coefficients
for calculating aqueous activity coefficients using the Truesdell-Jones equation.
The variables for all species are stored as the Matlab array called �databaseD�
which is called by the �fugCoefAQ� function described above. Note, however, that
this approach has been superseded by the Pitzer equation method.

load_Pitzer:
This function opens the pitzer_coef.txt file that contains parameters for
calculating aqueous activity coefficients using a simplified version of the
Pitzer equations. The variables for all species are stored in the Matlab
arrays B0, B1, B2, and C0 (binary interaction parameters), which are
called by function Pitzer_activity_diseq. This module goes unused for gas
phase calculations.
#~*~*~**~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~*~*

EXPLANATION OF INPUT FILES:

Two examples of input files are provided: inputs_Mars.txt and
inputs_ModernEarth_gas_only.txt
They both have the following structure:

- First line: Vector containing charge of all species. All gaseous species must
have zero charge.
- Second line: Vector containing phase identifiers for each species.
% Identifiers:
      0 - solvent (liquid water must be the first species in all vectors)
      1 - gases
      2 - liquids (databases for other liquids not included, but this is
included for a possible future expansion of the code)
      3 - solids (databases for solid species not included,
but this is included for a possible future expansion of the code)
      4 - aqueous electrolytes
- Third line: List of molecular species.
- Fourth line: Vector containing molar abundances for each species.
The molar abundances of gaseous species should be mixing ratios to
ensure the final Gibbs free energy change has units of joules per mole of
atmosphere. That is, the molar abundances of gaseous species should sum to 1.
- Fifth line: Vector containing the scalars ocean volume (om),
salinity scalar (ams), and a scaling factor (scale_factor).
The ocean volume can be changed to perform sensitivity tests using different
ocean volumes (but note that initial carbon chemistry must be re-equilibrated
for a fair comparison). The salinity scalar is no longer used. The scaling
factor helps ensure more rapid convergence of optimization algorithms.
It should not need to be changed for any of the examples provided. The scale
factor may not be necessary for python optimization but we include it in this
version of the model for consistency.
#~*~*~**~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~*~*

DATABASE EXPLANATION

NEWNASA:
The NEWNASA database contains the coefficients required to calculated the Gibbs
free energy of formation as a function of temperature for a large number of
molecular species. Appendix A in Krissansen-Totton et al. (2016) describes
how these calculations are performed. The format of the data entries in the NASA
text file is described here: http://www.grc.nasa.gov/WWW/CEAWeb/def_formats.htm.
The NASA database itself can be accessed here: http://www.grc.nasa.gov/WWW/CEAWeb/
whilst an updated version of the text file is available here:
"http://garfield.chem.elte.hu/Burcat/NEWNASA.TXT"


fugacityCoefficientsVariables:
This text file contains critical temperatures, critical pressures, and acentric
factors for a large number of gaseous species. There are used to calculate
fugacity coefficients using the Soave equation as described in Appendix A of
Krissansen-Totton et al. 2016. These thermodynamic parameters were obtained
from Knovel�s online database (see Perry, R. H., D. W. Green, and Knovel
(Firm) (2008), Perry's chemical engineers' handbook, 8th ed., 1 v. (various pagings)
pp., McGraw-Hill, New York.)


AQ_coefficients.txt:
This text file contains coefficients for calculating aqueous activity coefficients
using the Truesdell-Jones equation, as described in Appendix C of
Krissansen-Totton et al. (2016). The coefficients were sourced from Langmuir (1997).
However, this methodology has been superseded by calculating aqueous activity
coefficients using the Pitzer equations, which is a more accurate approach.


pitzer_coef.txt:
This text file contains the coefficients used for calculating aqueous activity
coefficients using the a simplified version of the Pitzer equations, namely the
binary interaction parameters B0, B1, B2, and C0. These interaction parameters were
obtained from Appelo and Postma (2005) and Marion (2002). See Appendix C in
Krissansen-Totton et al. (2016) and the Materials and Methods section in
Krissansen-Totton et al. (2018) for a detailed description of how these
parameters are used to calculate water activity and aqueous activity
coefficients, respectively.


sprons96_edited2.dat:
This file contains the sprongs96 SUPCRT database for calculating aqueous
Gibbs energies of formation. The aqueous species are clearly demarcated under the
header �aqueous species�. See Appendix C in Krissansen-Totton et al. (2016) for
equations describing how Gibbs energies of formation are calculated using the
species-specific coefficients provided in this database.
#~*~*~**~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~**~*~*~*~*~*

Contact email for code questions specific to the Python version:
amber.v.young@nasa.gov (preferred)
or
ambervyoung86@gmail.com

For questions regarding the Matlab version or general correspondece feel free
to contact either myself or Josh Krissansen-Totton (joshkt@uw.edu).
