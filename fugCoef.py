#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:51:18 2019

@author: avbritt
"""

"""
This module returns the value of the fugacity coefficient at the given
temperature and pressure for the set of gaseous species with the specified
mole amounts. This uses equations 27 and 28 in Krissansen-Totton et al. (2016)
Throws errors if the designated species cannot be located in the database or
if the computer has trouble solving the cubic equation for the compressability
factor.
"""

import numpy as np
import pandas

def fugCoef(temperature,pressure,names,n):
    #*****************************Constants***********************************
    #Ideal gas Constant (J/(mol*K))
    R = 8.314472
   #**************************************************************************
    #Load DatabaseB
    F = pandas.read_csv('fugacityCoefficientVariables1.txt', header=None,\
                        delimiter='\t', encoding='UTF-16BE')
    #Defining the column headers of the file. We have Species name, The
    #   Critical Temperature (in Kelvin), The Critical Pressure (in Mega
    #   Pascals), and finally the acentric value
    F.columns = ['Species','T_Crit','P_Crit','a']
    #Making the row index name the Species name, so the data for each species
    #   can be recovered using the species name as a key
    F = F.set_index('Species')
    #Multiplying the Critical Pressure column by 10 to convert from MPa to bar
    F['P_Crit'] = F['P_Crit']*10
    #Initialize empty lists for ai, bi, alphai, k
    ai = []
    bi = []
    alphai = []
    kij = 0
    #Use a for loop to loop through each species name
    for i in range(len(names)):
        #if names in F['Species']:
    #This takes a constant and multiplies it by the Critical Temperature of
    #   that species divided by the Critical Pressure of that species
        ai_i = 0.42747*(R**2)*F.loc[names[i]]['T_Crit']/F.loc[names[i]]\
        ['P_Crit']
        ai.append(ai_i)
        bi_i = 0.08664*R*F.loc[names[i]]['T_Crit']/F.loc[names[i]]\
        ['P_Crit']
        bi.append(bi_i)
        #
        alphai_i = (1+(0.48508+1.55171*F.loc[names[i]]['a']-0.15613*(F.loc[names[i]]['a']**2))*(1-np.sqrt(temperature/F.loc[names[i]]['T_Crit'])))**2
        alphai.append(alphai_i)

    #Define the variable representing the coefficient 'a alpha, total'.
    aaTotal = 0.
    aai = []
    #Over the length of a and alpha ... twice
    #Here were are computing the mixing rule of the Soave equation of state
    #   according to equation 28 of Krissansen-Totton et al. (2016)
    #   In the above paper they performed a sensitivity test for kij at temps
    #   and pressures high enough for a non-ideal gas regime. It turns out
    #   kij has little effect on the fugacities and therefore little effect
    #   on the overall Gibbs free energy. Therefore kij is taken to be 0 for
    #   every pair of species i and j. (In this case q=i and j=q+1)
    for q in range(len(ai)):
        for p in range(len(ai)):
            aai_i = (1-kij)*np.sqrt(ai[q]*alphai[q]*ai[p]*alphai[p])
            aai.append(aai_i)
            aaTotal = aaTotal + n[q]*n[p]/np.sum(n)**2*aai[q]

    #Calculate the coefficient representing 'a alpha, w/ indices'
        """if q == len(ai)-1:
            continue
        else:
            aai_i = (1-kij)*np.sqrt(ai[q]*alphai[q]*ai[q+1]*alphai[q+1])
            aai.append(aai_i)
    #Total the values to find the value of 'a alpha, total'
            aaTotal = aaTotal + n[q]*n[q+1]/np.sum(n)**2*aai[q]"""
    #Calculate the value of bTotal, the coefficient representing 'b w/o
    #   indicies'
    bTotal = []
    A = []
    BTotal = []
    Bi = []
    for i in range(len(n)):
        bTotal_i = (n[i]/sum(n))*bi[i]
        bTotal.append(bTotal_i)
    #Calculate the value of the coefficient representing 'A'
    #for i in np.arange(len(aaTotal)):
    for i in range(len(ai)):
        Ai = aaTotal*pressure/((R*temperature)**2)
        A.append(Ai)
    #Calculate the value of the coefficient representing 'B w/o indicies'
        BTotali = bTotal[i]*pressure/(R*temperature)
        BTotal.append(BTotali)
    #Calculate the value of the coefficient representing 'B w/ indicies'
        Bii = bi[i]*pressure/(R*temperature)
        Bi.append(Bii)
    #Coefficient terms to solve for compressability factor, z
    cubic = []
    for i in range(len(BTotal)):
        cubici = [1, -1, (A[i]-BTotal[i]-BTotal[i]**2), -A[i]*BTotal[i]]
        cubic.append(cubici)
        #print(cubic)
    #Use the customized method 'solveCubic' to calculate the value of z,
    #   the compressability factor
        zi = np.roots(cubic[i])
    #Function that returns all three roots of a function (The real root still
    #   needs to be solved for)

    #Initialize the compressability factor to infinity to help solve for the
    #   real root value
        z = np.inf
    #For loop for each root of z
    #Find the real root that is positive (z must be positive by definition)
    #   and is the closest value to 1.
        for q in range(3):
            if np.isreal(zi[q]) and zi[q]>0 and np.abs(zi[q]-1)<np.abs(z-1):
                #z = []
                #z_i = zi[q]
                #z.append(z_i)
                z = zi[q]
        if z == np.inf:
        #If no values of z fit that description then something is wrong
            print('Catling:badIncompressability',\
                  'The value of z is either non-real or less than zero.')

    #All coefficients have been found. Solve for the ln(Phi)
    #These are just smaller parts of the whole equation (27) in Krissansen-
    #   Totton et al. (2016)
    lnPhi = []
    for i in range(len(bTotal)):
        PhiA = Bi[i]/BTotal[i]*(z-1)
        PhiB = np.log(z-BTotal[i])
        PhiC = A[i]/BTotal[i]*(Bi[i]/BTotal[i]-2/aaTotal*(aai[i]*(n[i]/sum(n))))
        PhiD = np.log(1+BTotal[i]/z)
        lnPhii = PhiA-PhiB+PhiC*PhiD
        #AVY- Testing, lnPhi is assumed to be 0 here because it is neglibile
        #   at these temp-press, will add this back in once debugging is
        #   complete.
        lnPhi.append(lnPhii*0)
        """if len(z) == len(bTotal):
            PhiA = Bi[i]/BTotal[i]*(z[i]-1)
            PhiB = np.log(z[i]-BTotal[i])
            PhiC = A[i]/BTotal[i]*(Bi[i]/BTotal[i]-2/aaTotal[i]*(aai[i]*(n[i]/sum(n))))
            PhiD = np.log(1+BTotal[i]/z[i])
            lnPhii = PhiA-PhiB+PhiC*PhiD
            lnPhi.append(lnPhii)
        else:
            print('We may be missing compressability values')
            print(z)
            print(bTotal)
            print(Bi)"""

    Phi = [np.exp(i) for i in lnPhi]
    # print(lnPhi)

    #AVY- FOR Testing
    #return F, ai, bi, alphai, aai, aaTotal, lnPhi
    return lnPhi
