#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 18:52:58 2019

@author: avbritt
"""

"""
This module returns the value of the Gibbs free energy of formation for
aqueous species at the specified temperature and pressure. Throws errors if
the designated species cannot be located in the database or if the species
do not have an acceptable temperature range

This module uses equation 32 in Krissansen-Totton et al. (2016) to calculate
Gibbs free energies of formation. Parameters are taken from the sprons96
database.
"""

import numpy as np
import loadDatabaseC as ldC
import dielectric as de

def gibbsAQ(species1,T,P):
    #Load the database
    All_Data = ldC.loadDatabaseC()
    #Make it so that any terms (like +) translate literally when searching

    #Get the coefficient data for the species from the database
    if species1 in All_Data:
        coef_str = All_Data[species1]
    #Converting the list of sting type coefficients to floats
        coef = np.array(coef_str,dtype=float)

    #See equation 32 in Krissansen-Totton et al. (2016) for a definition of
    #   all the variables and coefficients below:
        Tr = 298.15
        Pr = 1
        Psi = 2600
        Theta = 228
        Y = -5.81*10**-5
        Gr = coef[0]
        Hr = coef[1]
        Sr = coef[2]
        a1 = coef[3]/10
        a2 = coef[4]*100
        a3 = coef[5]
        a4 = coef[6]*10000
        c1 = coef[7]
        c2 = coef[8]*10000
        w = coef[9]*100000
        q = coef[10]
        #Make a call to the dielectric module
        diE = de.dielectric(T,P)

        #These formulas are all parts of the formulas used in the Walther's
        #   Essentials of Geochemistry.
        G1 = Gr
        G2 = -1*Sr*(T-Tr)
        #np.log is the natural log
        G3 = -1*c1*(T*np.log(T/Tr)-T+Tr)
        G4 = a1*(P-Pr)

        h1 =np.log((Psi+P)/(Psi+Pr))

        G5 = a2*h1

        h2 = (1/(T-Theta))-(1/(Tr-Theta))
        h3 = (Theta-T)/Theta
        h4 = np.log(Tr*(T-Theta)/(T*(Tr-Theta)))

        G6 = -1*c2*(h2*h3-T/(Theta*Theta)*h4)
        G7 = (1/(T-Theta))*(a3*(P-Pr)+a4*h1)
        G8 = w*Y*(T-Tr)

        G = 4.184*(G1+G2+G3+G4+G5+G6+G7+G8)


    #If the species isn't there, throw a print statement
    #   NOTE: since the Database is saved as a dictionary, all key names must
    #   be unique and therefore, there can be no duplicate matches to
    #   species so no need to check for that
    else:
        print('Species',species1,'is missing from the aqueous database')
        G = ['G has not been found']


    return G
