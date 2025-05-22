#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:34:40 2020

@author: avbritt
"""

"""
This module returns the value of the Gibbs Free Energy of formation at 1 atm
for the given aqueous species at the specified temperature. Throws errors
if the designated species cannot be located in the database or the species
does not have an acceptable temeperature range.
"""

import numpy as np
import pandas as pd



import loadDatabase as ld
import makeAs as mAs
import convert2standard as c2s
import Parsing as Pars

#This shell function makes it so that the function can accept either a
#   cell array of multiple species, or a single species.
def gibbs(species,T):
     #*~*~*~*~*~*~*~*~*~* FUNCTIONS SECTION *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

    #gibbsHelper is a function to calculate the aparent Gibbs free energy
    #   of each species
    def gibbsHelper(species, T):
        #load the database
        #Note that NewNASA.txt database will return the Species coefficient
        #   data and the name dictionary since some species had the same
        #   name but different chemical structure (for these cases the data
        #   dictionary has arbitrary indexes to keep the species key
        #   names unique).
        All_Data, Data_names = ld.loadDatabase()

        #We need an if statement for the species NH3 in the Earth templates
        #   because NH3 has two entries with the same name in the NewNASA
        #   Database. The Matlab version of this code was using the second
        #   NH3 entry sourced from Gurvich, 1989 which is called NH3_2
        #   in this Database and is used by default.
        if species == 'NH3':
            species = 'NH3_2'
        #Find the matching species:
        if species in All_Data:
        #Goal is to loop over the columns of each species, taking the min
        #   assumes that there will inherently be less columns than rows
        #   which may or may not be the right approach but works for the time
        #   being
            a = np.min(np.shape(All_Data[species]))
        #initializing the coef array
            coef = []

        #A for loop to loop through the columns of the database
            for i in range(a):
        #If the given temperature is within the temperature range of the
        #   database (Temperature information is in the first two rows
        #   of the species data)
                # if T in All_Data[species][0:2,i]:
                if T > All_Data[species][0,i] and T < All_Data[species][1,i]:
        #Find the correct temperature range ans store that set of coefficients
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
            print('Species',species,'is missing from the gibbs gas database')
            G = ['G has not been found']
        return G
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~~*~*
    #Call to makeAs module which calculates the coefficient matrix 'a'
    #   (ie: CH4 is C1 H 4)
    #AVY NOTE: As of 4/17/2020 the math perfomred here to compute the Gibbs
    #   free energy of formation is computed one species at a time.
    a = mAs.makeAs(species)
    #print(a)


    #for loop that goes through each species
    #for i in range(len(species)):
    #Calculate the apparent gibbs free energy of that species
    #   (call gibbs and subfunction gibbsHelper which is already defined
    #   in the gibbsBB module)
    #If it's only one species, just calculate the gibbs free energy of the
    #   specific species
    G = gibbsHelper(species,T)
    #print(species,G)

    #Transpose the matrix so that the matrix algebra works later.
    G = np.array(G).T
    #Backup dataframe of the coefficient matrix in the form of a pandas
    #   datframe
    a1=pd.DataFrame.from_dict(a)
    #For Converting Berman-Brown Convention to Standard Gibbs free energies
    #  of formation
    #Call to the module convert2standard where we convert the dictionary that
    #   was created with makeAs into a new dataframe with each species in
    #   its standard state (if applicable, otherwise the species remains
    #   unchanged in the new dictionary)
    elementSpecies = c2s.convert2standard(a,T)
    #print(elementSpecies)

    #Create a cell array containing the number of each reference state species
    #   in each species

    #The gibbs free energy contribution by the elements for each elemental
    #   species
    element_g = []
    #Loop over each standard state element within the dictionary
    #   and compute the Gibbs free energy
    #   contribution from that element
    for element in elementSpecies:
    #Adjust that value by the number of moles each element has in it (H2 has
    #   2, for instance and the calculation below would reflect Gibbs_H2/2)
        G_element = gibbsHelper(element,T)/elementSpecies[element]
        # print(element,G_element)
        # print(element, G_element2)

    #Now we put the Gibbs energy of each element into the element_g array
        element_g.append(G_element)


    #With all those arrays, multiply those values by the number of each
    #   element the species have in them.
    #   (So, since CH4 has 4 H's, subtract
    #   the gibbs free energy of formation value by 4*(Gibbs_H2)/2)
    #Make sure element_g is in array form
    element_g = np.array(element_g)
    #@ is the shorthand for python dot product
    G = G - a1 @ element_g
    # print(G)

    return G
