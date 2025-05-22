#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:50:32 2019

@author: avbritt
"""

"""
This module returns the value of the acivity coeffiencient at the given
temperature and pressure for the set of given aqueous species with the
specified mole amounts. The Truesdell-Jones equation (described in
Krissansen-Totton et al. (2016) EQ. 33) is used to calculate aqueous activity
coefficients. This function has been superceded by the Pitzer equation method,
which is more accurate for Earth-like systems.

"""

import numpy as np
import loadDatabaseD as ld


def fugCoefAQ(temperature,pressure,names,n,v,om,scale_factor):
    #Constants in the Truesdell-Jones equation
    A = 0.5092 #Valid at 25 degrees C
    B = 0.3283 #Valid at 25 degrees C
    #loading the database with aqueous species information
    All_data = ld.loadDatabaseD()
    #Mass of the Ocean in kg
    mass_ocean = om*(1.34e9*1e9*1030)
    #Total moles of the atmosphere
    moles_atm = 1.7620e20
    #convert normalized moles back to molalities
    n_molal = [i*moles_atm/mass_ocean for i in n]
    # for i in range(len(n)):
    #     n_i = n[i]*moles_atm/mass_ocean
    #     n_molal.append(n_i)


    #I represents the ionic strength of the solution
    """for i in np.arange(len(n_molal)):
        I = 0.5*np.sum(n_molal[i]*v[i]*(v[i]/scale_factor))"""


    #Arrays used to hold the species specific AQ coefficients obtained from
    #   Langmuir 1997 (ie the info from AQ_coefficients1.txt)
    a_array = []
    b_array = []
    #Initializing new arrays that will only contain the charged, aqueous
    #   species
    names_aq = []
    n_aq = []
    v_aq = []

    #Loop through the charge vector and find all the species that have a
    #   non-zero charge
    for i in range(len(v)):
        if abs(v[i]) > 0:
    #If the species are in the same order as the charge vector, then the index
    #   of non-zero charge is at the same index of name of the species
            charged_species = names[i]
            n_aqi = n[i]
            v_aqi = v[i]
            names_aq.append(charged_species)
            n_aq.append(n_aqi)
            v_aq.append(v_aqi)

    #I represents the ionic strength of the solution
    I = 0.5*np.sum(n_molal*v*v/scale_factor)
    lnPhiAQ = np.zeros(len(names))
    #Loop through the names vector and find all the species that have a
    #   non-zero charge
    for i in range(len(names)):
    #If the species are in the same order as the charge vector, then the index
    #   of non-zero charge is at the same index of name of the species
        if v[i] != 0:
            lnPhiAQ[i] = -A*v[i]**2*(np.sqrt(I)/(1+B*All_data.loc[names[i]]['a']*np.sqrt(I)))+All_data.loc[names[i]]['b']*I
            # lnPhiAQ.append(lnPhiAQi)
        else:
    #lnPhiAQ is 0. if the species has no charge
            lnPhiAQ[i] = 0.

    return lnPhiAQ
