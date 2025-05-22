#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 18:40:30 2019

@author: avbritt
"""

"""
The function Gibbs_energy_minimization takes a (random) initial state
and calculates the multiphase equilibrium appyling atom and charge
conservation constraints, which are derived from the true initial state
(in load_vectors). It returns equilibrium abundances
and the Gibbs energy difference between the (true) initial
and equilibrium states. The only parameters that should be modified in
this file are temperature and pressure, and potentially some of the minimization
parameters used for the scipy minimize function if the optimization routine is
having trouble finding minimia.
"""

import numpy as np
import scipy.optimize as sciopt
from scipy.optimize import NonlinearConstraint, Bounds
import pandas as pd


import fugCoef
import fugCoefAQ
import load_vectors
import gibbsAQ as gAQ
import gibbs
import gibbsBB as gBB
import makeAs as mAs
import Parsing as Par
import Pitzer_activity_diseq as Pad


def Gibbs_energy_minimization(n):
    #Declaring "universial" parameters to be used throughout the module
    global R,T,P,l,v,g,V,names,n_true,a,n_true_preserved,om,scale_factor,b,c,\
        All_elements_dataframe,V,v_dataframe


    #Universal Gas Constant (J/K/mol)
    R = 8.3145
    #~*~*~*~*~*~*~*~*~*~*~*~*~INPUT REQUIRED*~*~*~*~*~*~**~*~*~**~*~*
    #Temperature of the system in Kelvin
    T = 288.15
    #AVY- For retrievals read in the temperature from file
    # T = np.genfromtxt('Retrieval_Temperature.txt',dtype=float)
    #Pressure of the system in Bars
    P = 1.013
    #AVY- For retrievals read in the pressure from file
    # P = np.genfromtxt('Retrieval_Pressure.txt',dtype=float)
    #~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

    #Load the input files from the load_vectors module
    #   name_input is an array of strings containing the names of each species
    #   l = array containing the phase of each species
    #   v = array containing the charge of each species
    v, l, n_init, om, sms, scale_factor, names = \
    load_vectors.load_input()


    #AVY-So currently the makeAs function works by returning a single species
    #   dictionary comprised of the constituent elements.
    All_elements_dict = {}
    for species in names:
        elements = mAs.makeAs(species)
        for dict_element in elements:
            if dict_element in All_elements_dict:
                All_elements_dict[dict_element].extend(elements[dict_element])
            else:
                All_elements_dict[dict_element] = elements[dict_element]
    # print(All_elements_dict)
    #Now that we have a dictionary of all the elements we can loop through
    #   the species list and elements in order to fill the dictionary such
    #   that it resembles a square matrix.

    #Dummy counter variable that is in sync with the species for loop
    #   to ensure the 0's in the dictionary are placed in the right spot.
    counter = 0
    #Dummy counter variable to skip over elements contained within the species
    #   since those species already have data entries in the dictionary
    element_counter = 0
    for species in names:
    #Call to the Parser module which returns a variable containing the
    #   list of elements for a particular species
        species_elements = Par.parser(species)
        # print(species, species_elements)
    #A for loop to loop over all elements within the element dictionary
        for dict_element in All_elements_dict:
    #A for loop to go over each element in the list that was returned
    #   from parser
            for single_element in species_elements:
                # print(dict_element,species_elements)
    #If the element is already in the list of elements for a particular
    #   species, skip over it by incrementing the dummy counter
                if dict_element in species_elements:
                    element_counter += 1
    #If the element is not within the list of elements for a particular
    #   species, add in a zero in the element dictionary.
                elif dict_element not in species_elements:
    #This if statement ensures we only insert a 0 once per not found
    #   species.
                    if single_element == species_elements[0]:
                        All_elements_dict[dict_element].insert(counter,0)
                  # print(dict_element,All_elements_dict[dict_element])
        counter += 1
    # print(All_elements_dict)

    #Now we can convert this dictionary (which now resembles a square matrix
    #   of element names acrocss the columns and the number of that
    #   element in each species along the rows)
    #   and convert it to a pandas dataframe.
    All_elements_dataframe = pd.DataFrame.from_dict(All_elements_dict)
    # All_elements_dataframe = All_elements_dataframe.set_index('names')
    #print(All_elements_dataframe)


    #~*~*~**~*~**~*~**~ FUNCTIONS SECTION ~*~**~*~*~*~**~**~**~*~*~*~*

    #This Function (onlyPhase) returns a form of some original matrix (old)
    #   with onlyspecies in it of the specified phase (phaseNumber)
    #Takes a given matrix and removes values from that matrix that are not
    #   assoiciated with the specified phase
    def onlyPhase(phaseNumber, old):
        global l
        new = []
        #For each species, if the species phase is the correct phase, add it to
        #   the new matrix (thus the new matrix will only contain species of
        #   the correct phase)
        for temp in range(len(old)):
            if l[temp] == phaseNumber:
        #The new version of the matrix with only the given phase
                new = np.append(new,old[temp])
        return new


    #The function totalGibbsInternal calculates the total Gibbs energy of the
    #   system given abundances, n, and the temperature and pressure of the
    #   system.
    #load in variables from the main function
    #   (R, T, g, P, l, names, v, n_true, scale_factor, om)
    #def totalGibbsInternal(R,T,g,P,l,names,v,n_true,scale_factor,om,n):
    def totalGibbsInternal(n):

        #Calculate the total moles of gases, liquids and solids by calling
        #   onlyPhase function which is defined in this module
        TotalMolesGases = np.sum(onlyPhase(1,n))
        TotalMolesLiquids = np.sum(onlyPhase(0,n))+np.sum(onlyPhase(2,n)) \
            + np.sum(onlyPhase(4,n))
        TotalMolesSolids = np.sum(onlyPhase(3,n))

        #The natural log of the fugacity coefficient for each gas phase species
        #   in the system
        lnPhi = fugCoef.fugCoef(T,P,onlyPhase(1,names),onlyPhase(1,n))
        lnPhi = np.array(lnPhi)
        # if len(onlyPhase(4,n)) != 0:
        #     lnPhiAQ = fugCoefAQ.fugCoefAQ(T,P,onlyPhase(4,names),onlyPhase(4, n),\
        #                         onlyPhase(4,v),om,scale_factor)
        #     lnPhiAQ = np.array(lnPhiAQ)

        #fugCoefAQ outputs are superceded by Pitzer_activity_diseq outputs
        #   (see below), but the function is retained in case the user wishes
        #   explore regimes where the Truesdell-Jones equation is more
        #   to appropriate than Pitzer equations.


        #Calculate water activity and aqueous species activity coefficients
        #   using Pitzer equations
        #Function call to Pitzer_activity_diseq() which is a seperate module
        # lnPhiWater = Pad.Pitzer_activity_diseq(names,n,v,l)
        #Fill in activity coefficients for anions

        #Fill in activity coefficients for cations

        #AVY- Need to properly calculate the fugacity of water using the
        #   pitzer equation module but for now just use a dummy value for
        #   testing Using a sample value calculated in the Matlab version for
        #   the Archean max template
        # lnPhiWater =   -0.019791506335140

        #The contribution to the total Gibbs free energy from all the species
        #   of each phase of the system.
        #This is essentially equation 10 from Krissansen-Totten et al. (2018):
        funct = np.zeros(2)

        #List with the abundaces of only the solvent species
        #   Inherent assumption is that the solvent is water since we use
        #   the pitzer_activity module to calculate the activity of
        #   water specifically (and additionally, the aqueous species as well)
        # solvent = onlyPhase(0, n)

        #If the solvent phase is not empty perform the calculation below
        # if len(solvent) != 0:
        #     funct[0] = onlyPhase(0,n).T*(onlyPhase(0,g)/R/T+lnPhiWater)
        # #print(funct[0])
        # #Otherwise the solvent Gibbs contribution is 0
        # else:
        #     funct[0] = 0.0

        #List with the abundaces of only the gas phase species
        gas_phase = onlyPhase(1, n)

        #If the gas phase is not empty perform the calculation below
        if len(gas_phase) != 0:
        #For loop to calculate the Gibbs free energy contribution of the
        #   gas phase speices at each fugacity (lnPhi) value
            # for i in range(len(lnPhi)):
                # funct_val = np.sum(onlyPhase(1,n).T*(onlyPhase(1,g)/R/T+np.log(P)+\
                #                         lnPhi[i]+np.log(onlyPhase(1,n)/TotalMolesGases)))
            funct_val = np.sum(onlyPhase(1,n)*(onlyPhase(1,g)/R/T+np.log(P)+\
                                    lnPhi+np.log(onlyPhase(1,n)/TotalMolesGases)))

            #print(funct_val)
            funct[1] = funct[1] + funct_val
            #print(funct[1])
        #Otherwise the gas Gibbs contribution is 0
        else:
            funct[1] = 0.0
        #List with the abundaces of only the liquid phase species
        liquid_phase = onlyPhase(2,n)

        #If the liquid_phase list is not empty then perform the the calculation
        #   for the Gibbs free energy contribution from the liquid phase species
        # if len(liquid_phase) != 0:
        #     funct[2] = np.sum(onlyPhase(2,n)*(onlyPhase(2,g)/R/T+np.log(P)+\
        #                             np.log(onlyPhase(2,n)/TotalMolesLiquids)))
        # #Otherwise the contribution from that phase is 0.
        # else:
        #     funct[2] = 0.0
        # #print(funct[2])
        #List with the abundaces of only the solid phase species
        # solid_phase = onlyPhase(3, n)
        #If the aolid_phase list is not empty then perform the the calculation
        #   for the Gibbs free energy contribution from the solid phase species
        # if len(solid_phase) != 0:
        #     funct[3] = np.sum(onlyPhase(3,n).T*(onlyPhase(3,g)/R/T+np.log(P)+\
        #                             np.log(onlyPhase(3,n)/TotalMolesSolids)))
        # #Otherwise the contribution from that phase is 0.
        # else:
        #     funct[3] = 0.0

        #norm_liq = onlyPhase(4,n)/TotalMolesLiquids
        #log_norm_liq = [np.log(i) for i in norm_liq]
        #print(log_norm_liq)
        #List with the abundances of only the aqueous phase species
        # aqueous_phase = onlyPhase(4, n)
        # if len(aqueous_phase) != 0:
        #For loop to calculate the Gibbs free energy contribution of the
        #   aqueous phase speices at each fugacity (lnPhiAQ) value

        # funct_val = np.sum(onlyPhase(1,n)*(onlyPhase(1,g)/R/T+np.log(P)+\
        #                             lnPhi+np.log(onlyPhase(1,n)/TotalMolesGases)))

            # for i in range(len(lnPhiAQ)):
                #print([np.log(i) for i in onlyPhase(4,n)/TotalMolesLiquids])

            # funct_val_aq = np.sum(onlyPhase(4, n)*(onlyPhase(4,g)/R/T+lnPhiAQ+np.log(55.508435)+\
            #                       np.log(onlyPhase(4,n)/TotalMolesLiquids)-\
            #                           np.log(onlyPhase(0,n)/TotalMolesLiquids)))
            # funct[4] = funct[4] + funct_val_aq
        #Otherwise the contribution from that phase is 0.
        # else:
        #     funct[4] = 0.0
        #Here the contribution by each phase are totaled together into a total
        #   Gibbs free energy of the system
            #print(funct)
        fun = 0

        #Calculating the total Gibbs free energy contributions of each species
        # as a function of phase
        fun = np.sum(funct)

        #Calculate analytic gradient of the system for each species in the
        #  system: This follows equation 36 in Krissansen-Totten et al. (2016)
        # AVY - This particular scipy minimizer does not require the user to
        #   compute the gradiant but the code snipets below are left in as a guide for
        #   including it if the user wishes to do so or if they choose a different
        #   minimize routine to use.
        # grad = np.zeros(len(n))
        # #For each species, identify its phase of matter and use the appropriate
        # #   formula
        # for i in range(len(n)):
        #     #for j in range(len(l)):
        #     if l[i] == 0:
        #         grad[i] = g[i]/R/T+np.log(n[i]/TotalMolesLiquids)-\
        #         TotalMolesLiquids/n[i] - n[i]/TotalMolesLiquids + 2.
        #     elif l[i] == 1:
        #         grad[i] = g[i]/R/T+np.log(P)+lnPhi[0]+np.log(n[i]/TotalMolesGases)
        # #Removes the first fugacity coefficient from the list of coefficients.
        # #   That way, each coefficient is only used once.
        #         lnPhi = lnPhi[1:]
        #     elif l[i] == 2:
        #         grad[i] = g[i]/R/T+np.log(P)+np.log(n[i]/TotalMolesLiquids)
        #     elif l[i] == 3:
        #         grad[i] = g[i]/R/T+np.log(P)+np.log(n[i]/TotalMolesSolids)
        #     elif l[i] == 4 and lnPhiAQ:
        #         grad[i] = g[i]/R/T + np.log(55.508435) + lnPhiAQ[0] + \
        #             np.log(n[i]/TotalMolesLiquids) - np.log(n[0]/TotalMolesLiquids) -\
        #             n[0]/TotalMolesLiquids + 1.
        #         lnPhiAQ = lnPhiAQ[1:]

        return fun

    #This Function (mycon) is used by the minimization function to represent
    #   the charge balance
    def mycon(x):
        #v = charges for each molecule; V = total charge of the system
        #   (should be 0 with sensible inputs)
        global v, V
        #Compute nonlinear inequalities at x. This is blank since we have no
        #   inequalities
        c = []
        #Compute nonlinear equalities at x (the moles values)
        ceq = (v*x-V)
        cd = []
        ceqd = np.array(v).T
        return c, ceq, cd, ceqd

    #This function is used by the minimization routine
    #   as part of the mass constraints we put on the system where
    #   A*n' = A*n.
    #   Here A is the dataframe of elemements, n' represents the new species
    #   mixing ratios and n represents the initial species mixing ratios
    #   Here c = A*n and b = A*n'
    def mass_constraint(n):
        global All_elements_dataframe, c

        b = All_elements_dataframe.T@(n)
        b = pd.DataFrame(b)
        # print(c)
        # print(b)
        # print((c - b).values.flatten())
        # print(n)

        return (c - b).values.flatten()


    #This function is used by the minimization routine to implement
    #   charge constraints on the system where
    #   v*n' = v*n
    #   Here v is the dataframe of the charge on each species, n' represents
    #   the new species mixing ratios and n represents the initial species
    #   mixing ratios
    #   Here V = v*n and VV = v*n'
    # def charge_constraint(n):
    # #   V = v_dataframe.T@(n_true_dataframe)
    #     global v_dataframe, V
    #
    #     VV = v_dataframe.T@(n)
    #     # print(VV)
    #     # print(v_dataframe)
    #     # print(V)
    #     VV = pd.DataFrame(VV)
    #     # print((V-VV))
    #
    #     return (V-VV).values.flatten()
    #Combined Contraints
    # def combined_const(n):
    #     global v_dataframe, V, All_elements_dataframe, c
    #     x = []
    #     VV = v_dataframe.T@(n)
    #     VV = pd.DataFrame(VV)
    #
    #     b = All_elements_dataframe.T@(n)
    #     b = pd.DataFrame(b)
    #     combined_c = np.append((c - b).values.flatten(),(V-VV).values.flatten())
    #     combined_c = np.abs(10**8-combined_c)
    #     # print(combined_c)
    #     return combined_c

    #This function is used by the minimization routine
    #   in order to specify constraints
    # def const(n):
    def const(n):
        #Implementing the mass constraint that A*n' = A*n
        #   where A is the dataframe of elements, n' represents the new species
        #   mixing ratios and n represents the initial species mixing ratios
        #   Here c = A*n and b = A*n'

        # if len(onlyPhase(4, n)) != 0:
        #     # constraint = [{'type':'ineq', 'fun': mass_constraint}, \
        #     #              {'type':'ineq', 'fun': charge_constraint}]
        #     constraint = {'type':'eq', 'fun': combined_const}
        #
        # else:
        constraint = {'type':'eq', 'fun': mass_constraint}
        # cons+=(constraint,)
        # constraint =
        # NonlinearConstraint(constraint,-np.inf,ub=1.9)

        return constraint
    # #Function to calculate the Jacobian
    # def jacob(n):
    #     global All_elements_dataframe, v

    #     x = np.append(v,All_elements_dataframe.T.values.flatten())
    #     return x

    # class atom_const:
    #     def _init_(self,n):
    #         global All_elements_dataframe, b
    #         self.const = b - All_elements_dataframe.T@(n)

    #~*~*~**~*~*~* END OF THE FUNCTIONS SECTION ~*~*~**~*~**~*~**~*~**~*~*~**~


    #Loop to calculate the Gibbs Free Eneergy of each species
    #AVY- Note that the functions gibbs, gibbsBB, and gibbsAQ
    #   are their own modules
    g = []
    for i in range(len(names)):
    #If the species are aqueous, use the aqueous database. Otherwise use the
    #   gaseous database. (NOTE: aqueous and multiphase calculations are not
    #   supported by the current version of the model.)
        if l[i] == 4:
    #Function Call to gibbsAQ which is a function that returns the aqueous
    #   Gibbs free energy values of formation
            g_val = gAQ.gibbsAQ(names[i],T,P)
            g = np.append(g,g_val)
        else:
    #Function call to gibbs and gibbsBB which returns the Gibbs
    #   free energy of formation for gaseous species
    #   This formulation is required to ensure the gaseous Gibbs energies of
    #   formation match the conventions for the aqueous species Gibbs energies
    #   of formation
            g_val = gibbs.gibbs(names[i],298.15)+gBB.gibbs(names[i],T)\
                - gBB.gibbs(names[i],298.15)
            g = np.append(g,g_val)

    #Use the true initial state to calculate atom and charge conservation
    #   conditions
    n_true = np.array(n_init).T

    #Option to force moles in atmosphere to sum to one (not recommended)


    #Convert molalities to moles per mole atmosphere for AQ species
    #Number of moles in ocean
    # mass_ocean = om*1.380200000000000e+21
    #
    # #Number of moles in atmosphere
    # moles_atm = 1.762035714285714e+20
    # #Find the indices for aqueous species
    # for i in range(len(names)):
    #     if l[i] == 4:
    # #Convert to moles per mole of atmosphere
    #         n_true[i] = n_true[i]*mass_ocean/moles_atm
    # #Ensure the net charge is 0 by adding Na(+)
    # Index_Na = np.where(np.asarray(names) == 'Na(+)')
    # #Adjust Na(+) to ensure zero net charge
    # n_true[Index_Na]=n_true[Index_Na]-v@n_true
    #
    # n_true_preserved = n_true

    #Multiply aundances by scaling factor for more precise calculations at low
    #   abundances
    n_true = [i* scale_factor for i in n_true]
    #Converting n_true to a Pandas Dataframe for matrix multiplication
    #All_elements_dataframe = pd.DataFrame.from_dict(All_elements_dict)
    n_true_dataframe = pd.DataFrame(n_true)


    #Charge conservation constraint (total charge should be close to zero)
    #   that will be used in the minimization routine
    v_dataframe = pd.DataFrame(v)
    # print(v_dataframe)
    # print(n_true_dataframe)
    V = v_dataframe.T@(n_true_dataframe)
    # print(V)

    #Calculate the total Gibbs energy of the true initial state
    #   Function call to totalGibbsInternal which is defined within this
    #   module and n_true is passed in as an arguement
    #gradient = []
    #G1, gradient = totalGibbsInternal(n_true)
    G1 = totalGibbsInternal(n_true)


    #atom conservation constraint that will be used in minimization routine
    c = All_elements_dataframe.T@(n_true_dataframe)
    # b = All_elements_dataframe.T@(n_true_dataframe)



    # #~*~*~**~*~**~*~*~* ADJUSTABLE VALUES ~*~*~**~*~*~***~*~*~*~*~*~*~**~*~*~
    # #   These values can be modified to change the accuracy of the algorithm
    # #AVY- May not be relevant for the Scipy function
    #
    # #The required tolerance in the Gibbs free energy before exiting the
    # #   minimization algorithm
    # tolfun = 1e-12 #1e-60
    # #The required tolerance in n values before exiting the minimization
    # #   algorithm
    # tolX = 1e-12 #1e-60
    # #Smallest possible real number for the input parameters
    # tolMin = 1e-12 #1e-60
    # #Maximum value of the constraint function
    # tolCon = 1e-12 #1e-60
    #
    # #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~*~*~*~*~*~*~*~*~*~


    #Run the minimization routine, starting from random initial state (n)
    #   Using atom conservation constraints, We're going to use
    #   scipy.optimize.minimize using the
    #   SLSQP method which uses least squares programming to minimize a function
    #   containing several variables. It includes bounds, equality and inequality
    #   constraints as specified by the user.
    #   We also use "options" to set the max number of iterations to 2000
    #   and we use "constraints" to pass our constraint dictionary to the
    #   solver.
    G2 = []
    outputs = []

    outputs = sciopt.minimize(totalGibbsInternal, n, method='SLSQP',\
                              options={'disp':True,'maxiter':2000},constraints=const(n),
                              bounds=Bounds(1e-308,np.inf,True))

    #The minimization routine will return an object,
    #   which is saved in the outputs variable.
    #   From outputs we obtain the equilibrium abundances (n),
    #   and total Gibbs Energy of the final state (G2),
    n = outputs['x']
    #Define G2 which is the total Gibbs free energy of the final state
    G2 = outputs['fun']


    #Define dG which is the change in Gibbs free energy of the system between
    #   initial and final state
    dG = G2 - G1

    #Calculation for the mass balance roughness
    #str(np.sqrt((b-All_elements_dataframe*n)*(b-All_elements_dataframe*n)))
    n_dataframe = pd.DataFrame(n)
    MBR = np.sqrt((c-All_elements_dataframe.T@(n_dataframe))\
                         *(c-All_elements_dataframe.T@(n_dataframe)))

    #Calculation for the Charge Balance
    #str(V-v*n)
    # CB = V-v_dataframe.T.dot(n_dataframe)
    # print(CB.iloc[0][0])

    #Output selected results.
    print('Mass Balance Roughness (should be zero): ',str(MBR.iloc[0][0]))
    # print('Charge Balance (should be zero): ',str(CB.iloc[0][0]))
    print('dG value: ',str(dG))
    print('G1 value: ',str(G1))
    print('G2 value: ',str(G2))
    print('deltaG value this iteration (J/mol): ',str(dG*8.3145*T/scale_factor))

    #Fill vals array with selected outputs
    vals = []
    for zz in range(len(names)):
        vals_value = n[zz]
        vals = np.append(vals,vals_value)


    vals = np.append(vals,dG)
    vals = np.append(vals,G2)


    #Display the initial and final abundances for every species
    #   AVY- Comment in for more detailed output
    # for ii in range(len(names)):
        # print(names[ii], '',n_true[ii]/scale_factor, '', n[ii]/scale_factor)
    return vals
