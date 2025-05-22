#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:59:09 2020

@author: avbritt
"""

"""
This module Pitzer_activity_diseq calculates the activity coefficients of
aqueous species and liquid water activity using a simplified version of the
Pitzer equations.

This module takes in the names of species, their abundances, charges, and
phases as inputs.
"""

import numpy as np
import load_Pitzer as lP

def Pitzer_activity_diseq(names_in,n_in,v_in,l):
    #Ensure Pitzer equation binary interaction parameters are loaded from
    #   the database file
    global B0, B1, B2, C0, om, scale_factor

    #Loading the Pitzer database
    #   Access data via dictionary keys B0, B1, ...
    All_data_Pitzer = lP.load_Pitzer()
    #Make n a row vector instead of a column vector (AVY-not sure if this
    #   is essential)
    n_in = np.array(n_in).T

    #Mass of the ocean in kg
    mass_ocean=om*1.34e9*1000**3*1030
    #Total moles of the atmosphere
    moles_atm=1.7620e+20
    #Convert the normalized moles back to molalities (meaning moles/kg)
    n=n_in*moles_atm/mass_ocean/scale_factor
    #Firstly, reduce names, abundance, and charge vectors to aqueous species
    #   only (select l=4 phase)
    #   Saving the first entry of the np.where tuple
    AQ_species_index = np.where(np.array(l) == 4)[0]
    names = []
    v = []
    m = []
    for index in AQ_species_index:
        names.append(names_in[index])
        v.append(v_in[index])
        m.append(n[index])

    #remove parentheses from species names to match the entries in the database
    for i in range(len(names)):
        names[i] = names[i].replace('(','').replace(')','')

    #Find the indices of all the positively charged (v>0) and
    #   negatively charged (v<0) aqueous phase species
    #   Saving the first entry of the np.where tuple
    cation_index = np.where(np.array(v)>0)[0]
    anion_index = np.where(np.array(v)<0)[0]

    #Create vectors containing cation names, charges, and abundances
    names_c = []
    v_c = []
    m_c = []
    for index in cation_index:
        names_c.append(names[index])
        v_c.append(v[index])
        m_c.append(m[index])

    #Create vectors containing the anion names, charges, and abundances
    names_a = []
    v_a = []
    m_a = []
    for index in anion_index:
        names_a.append(names[index])
        v_a.append(v[index])
        m_a.append(m[index])

    #The terms below are defined in equation 11 in Krissansen-Totton et al.
    #   (2018), originally from Marion and Kargel (2007)
    sum_mi = np.sum(m)
    #Ionic strength
    I=0.5*np.sum(np.array(m)*np.array(v)**2)
    #This could be modified to include temperature dependence, but is constant
    #   here (see Marion and Kargel (2007))
    A_phi=0.3915
    #Constant (see Marion and Kargel (2007))
    b=1.2
    Z=np.sum(np.array(m)*abs(np.array(v)))
    #Part of equation 11 in Krissansen-Totton et al. (2018)
    f_gamma = -A_phi*(I**0.5/(1+b*I**0.5)+2*np.log(1+b*I**0.5)/b)
    # print(f_gamma)

    #Create empty arrays for log activity coefficients
    ln_gamma_ct = np.zeros(len(v_c))
    ln_gamma_an = np.zeros(len(v_a))
    F=f_gamma
    #In the portions of the code that follow, B0, B1, B2, and C0 will be
    #   treated seperately. Note that B0, B1, and B2 are abbreviated forms
    #   of the terms B(0)_MX, B(1)_MX, and B(2)_MX, refferred to in KT16 and
    #   KT18. Bracketed numbers should be superscripted. C0 is a special case
    #   (see below)
    B0_sum = 0
    B1_sum = 0
    B2_sum = 0
    C0_sum = 0
    #Begin B0 Section
    #Find all cations in the first column of the database that match the
    #   species in our cation list.
    # ind1 = []
    # for i in range(len(names_c)):
    #     ind1.append(np.where(names_c[i]==All_data_Pitzer['B0'])[0])
    # print(ind1)

    #Find all anions in the second column of the database.

    #C contains indices of all cation-anion database rows relevant to binary
    #   interactions in this system.


    #Fill B0_used with relevant rows above: cation, anion, parameter, parameter
    #   parameter (last two columns will later be filled with abundances)

    #double summation in equation 34 in Krissansen-Totton et al. (2016)
    summation = B0_sum+B1_sum+B2_sum+C0_sum
    osmo = 1+(2./sum_mi)*(-A_phi*I**1.5/(1+b*I**0.5)+summation)

    #Activity of water
    a_w = np.exp(-osmo*sum_mi/55.50844)

    #take natural log before returning water activity
    lnPhiWATER=np.log(a_w)

    return lnPhiWATER
