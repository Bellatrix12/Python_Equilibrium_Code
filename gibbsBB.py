#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:30:02 2020

@author: avbritt
"""

"""
This module returns the value of the Gibbs Free Energy of Formation at 1 atm
for the given gaseous species at the specified temperature. Thows errors if
the designated species cannot be located in the database or the species does
not have an acceptable temperature range.

This module is similar to the gibbs.py module except that the Berman-Brown
convention for the Gibbs energies of formation is retained
"""

import numpy as np
import loadDatabase as ld

def gibbs(species,T):
    #*~*~*~*~*~*~*~*~*~* FUNCTIONS SECTION *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

    #gibbsHelper is a function to calculate the aparent Gibbs free energy
    #   of each species
    def gibbsHelper(species, T):
        #load the database
        #Note that NewNASA.txt database will return the Species coefficient
        #   data and the dictionary name since some species had the same
        #   name but different chemical structure (for these cases the data
        #   dictionary has arbitrary indexes to keep the species key
        #   names unique).
        All_Data, Data_names = ld.loadDatabase()
        #We need an if statement for the species NH3 in the Earth templates
        #   because NH3 has two entries with the same name in the NewNASA
        #   Database. The Matlab version of this code using the second
        #   NH3 entry sourced from Gurvich, 1989 which is called NH3_2
        #   in this Database and used by default.
        if species == 'NH3':
            species = 'NH3_2'
        #Find the matching species:
        if species in All_Data:
        #Goal is to loop over the columns of each species, taking the min
        #   assumes that there will inherently be less columns than rows
        #   which may or may not be the right approach but works for the time
        #   being
            a = np.min(np.shape(All_Data[species]))
        #initializing the coefficient array
            coef = []
        #A for loop to loop through the columns of the database
            for i in range(a):
        #If the given temperature is within the temperature range of the
        #   database (Temperature information is in the first two rows
        #   of the species data)
                if T > All_Data[species][0,i] and T < All_Data[species][1,i]:
        #Find the correct temperature range and store that set of coefficients
        #   into the variable 'coef'.
                    coef = All_Data[species][2:,i]
        #If none of the temperature ranges fit the designated temperature
            if len(coef) < 2:
        #If the temperature is less than the lowest temp
                if T < All_Data[species][0,0]:
        #Assign coef to be the lowest temperature range possible (which will
        #   hopefully be close).
                    coef = All_Data[species][2:,0]
                else:
        #Otherwise, use the highest range (since that should be closest to the
        #   higher temperature).
                    coef = All_Data[species][2:,-1]
        #Turn the coef array into an numpy array as opposed to a list
            coef = np.array(coef)
        #If the coefficients are in an acceptable form, calculate G.
            H = -1*coef[0]*T**(-2) + coef[1]*np.log(T)/T + coef[2] + coef[3]*T/2 + coef[4]*T**2/3 + coef[5]*T**3/4 + coef[6]*T**4/5 + coef[8]/T
            S = -coef[0] * T**-2/2 - coef[1]*T**-1 + coef[2]*np.log(T) + coef[3]*T + coef[4]*T**2/2 + coef[5]*T**3/3 + coef[6]*T**4/4 + coef[9]
            G = 8.3145*T * (H - S)
        #    print(G)
        else:
            print('Species',species,'is missing from the gibbsBB gas database')
            G = ['G has not been found']
        return G
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~~*~*


    #Call to makeAs module which calculates the coefficient matrix 'a'
    #   (ie: CH4 is C1 H 4)

    #Call to iscellstr() function which must be another Matlab specific
    #   function. If it's a list, then perform the function on them all.


    #for loop over the length of each species
    #G = []
    #for i in np.arange(len(species)):
    #Call to function gibbsHelper to calculate the apparent gibbs free energy
    #   of that species
    #   This retains the Berman-Brown convention
    G = gibbsHelper(species,T)

    return G
