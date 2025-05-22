#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 17:24:37 2019

@author: avbritt
"""

"""
This script calculates the available Gibbs free energy for the specified
single or multiphase system by repeatedly minimizing the total Gibbs energy of
the system starting from several (random) initial conditions.

Note: at the time of writing only the gas phase calculations can be run at this
time. Incorporating the multiphase calculations will be a product of future
versions of the Python model.
"""
#import Statements (built-in modules)
import numpy as np
#import random

#import Statements (our modules)
import load_vectors
import Gibbs_energy_minimization as Gem
import Plot_outputs as pltout

#Variable vals (empty array)
vals = []
#Set the number of iterations (num)
#AVY - For gas phase calculations only 1 iteration of num should be needed
#   to find the minimum Gibbs free energy value
num = 1

#A Call to load vectors
v_init, l_init, n_init, om, sms, scale_factor, name_input = \
    load_vectors.load_input()
#Multiply the initial abundances by a large scale factor for a more reliable
#   minimization
#   The next routine will divided the abundances by scale_factor at the end to
#   return to sensible units
n = [i * scale_factor for i in n_init]


#AVY- For coupling to retrievals, comment in these two lines below.
# vals_temp = Gem.Gibbs_energy_minimization(n)
# vals = np.append(vals,vals_temp)

#AVY - For Retrievals comment out this section *************************
#Code Segment for generating randomized initial conditions
#   Run the Gibbs free energy minimization routine and store the output in the
#   vals array (call to randperm)
for i in range(num):
    n = np.random.permutation(n)
    for j in np.arange(len(n)):
        #This generates totally arbitrary random initial conditions for each
        #   iteration
        n[j] = n[j]*500*np.random.rand()*10**(10*np.random.rand()-9)
    #Runs the Gibbs energy minimization routine and stores the output in the
    #   vals array
    vals_temp = Gem.Gibbs_energy_minimization(n.T)
    vals = np.append(vals,vals_temp)
#************************************************************************

#Find the global minimum value (using min function)
# Matlab version: [M,I] = min(vals[:,len(n)*2+3])
Min_val = min(vals)
#Comment in for more detailed output
# print(Min_val)

#Make an array of values divided by a scale factor
#Matlab version: array_of_Gibbs = vals[:,len(n)*2+3/scale_factor]
array_of_Gibbs = [(i * len(n)*2+3)/scale_factor for i in vals]
#Display (print) the Gibbs energy difference from all iterations
#Comment in for more detailed output
# print('Array containing available Gibbs energies from each iteration: ')
# print(array_of_Gibbs)

#We need to get the initial state vectors for plotting purposes, store the
#   initial values in new variables after calling (load vectors)
v_init, l_init, n_init, om, sms, scale_factor, name_input = \
    load_vectors.load_input()
n_true = n_init
l = l_init


#AVY - Aqueous sections go unused for gas phase calculations
#Change molalities to moles per mole of atmosphere for AQ species
#mass of ocean in kg
mass_ocean = om*1.3802e21
#total number of moles in Earth's atmosphere
moles_atm = 1.762035714285714e+20
#Find the indices for aqueous species only
aqueous_species = []
aqueous_index = []
for i in range(len(name_input)):
    if l[i] == 4:
        aq_species = name_input[i]
        index_val = i
        aqueous_species = np.append(aqueous_species,aq_species)
        aqueous_index = np.append(aqueous_index,index_val)
#Convert molalities to moles per mole atmosphere
for index in range(len(aqueous_index)):
    n_true[index] = n_true[index]*mass_ocean/moles_atm

#Extract final abundances for the global minimum case:
#   final_n is the final abundances for the global minimum case
min_val = len(n_true)+1
max_val = min_val+len(n_true)-1
#Matlab version: final_n = vals[I,min_val:max_val]/scale_factor
final_n = [i/scale_factor for i in vals]

#Find the index for the Gibbs Energy difference
G_dex = max_val+3
#   Print the largest (negative) Gibbs Energy difference
#Comment in print statements below for more detailed output
# print('Global minimum (J/mol):')
deltaG_value = vals[-2]/scale_factor
# print(deltaG_value)
#   Check if the Gibbs energy change is positive and print a warning if it is
#   global minimum has not been obtained, run the program again and try more
#   iterations
if deltaG_value > 0:
    print('Warning, positive Gibbs energy change - global minimum has not been obtained')
    print('Run the program again and try more iterations.')

#Call a function to plot the outputs (pass to the function n_true and final_n)
pltout.Plot_outputs(n_true,final_n)
