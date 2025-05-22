#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:18:08 2019

@author: avbritt
"""

"""
This module loads the database in the file 'F' into a dictionary.
"""
#Grabs whatever directory you are in (make sure you are in the top
#level directory of Multiphase_Equilibrium_Code)
import os
dirpathdefault = os.getcwd()

import numpy as np

def loadDatabase():
    #The name of the file holding the data
    #Note that a Unicode error showed up for the CH4+ species so we ignore
    #   the error in the open statement. Ignoring the error does not remove
    #   the species from the database (thankfully)
    with open('NEWNASA.txt','r+', errors='ignore') as F:
    #While loop for skipping introductory comments
        while F.readline().strip() != 'FOR 9 CONSTANTS SPECIES NOT INCLUDED IN THIS FILE SEE <http://cea.grc.nasa.gov>':
            continue
    #Skiping empty line
        _=F.readline()
    #Create and empty dictionary of species (Note that species with
    #   duplicate names but different chemical structures will be indexed
    #   as species_#. For example: Br2O_2)
        NEWNASA_Species = {}
    #Creating an empty dictionary that will only contain species names
    #   This is because some species have the same name but different
    #   chemical structure and its an easy way to reference what the index
    #   number means.
        NEWNASA_SPECIES_NAMES = {}
     #A variable to indicate when we've reached the end of the file,
     #In this case we reach the end of the file when 3 empty lines are found
        empty_count = 0
        while empty_count < 3:
    #line containing species information (whole line is a string)
            line_str = F.readline()
    #Splitting the line string into an array
            line = line_str.split()
    #If else statement to check to see if we are at the end of the file using
    #   the variable empty_count
            if line:
                empty_count = 0
            else:
                empty_count+=1
                continue
    #First character of the line
            First_Char =line[0].strip()[0]
            #if First_Char.isalpha():
    #If any of the characters within the line are letters, Then we're
    #   on the line that contains the species name information
            if any(c.isalpha() for c in line_str):
    #The species name is the first item in the line list
                species_name = line[0]
    #If statement to check if the species name line has phase information
                if line[1] in ['SOLID', 'LIQUID','GAS']:
    #Adding the phase to the species dictionary key name
                    species_name = species_name+'_'+line[1]
    #If statement to cheeck if there are duplicate species
                if species_name in NEWNASA_Species:
    #While loop for species that have multiple entries but
    #   different chemical structures. num is the index number being tacked
    #   on to the species name
                    num = 2
                    while f'{species_name}_{num}' in NEWNASA_Species:
                        num+=1
                    species_name = f'{species_name}_{num}'
    #Adding the line with the species name information to  it's own
    #   dictionary
                NEWNASA_SPECIES_NAMES[species_name] = line_str.strip()
    #Skipping the line after the species name because that information is not
    #   needed
                _=F.readline()
    #Check to see if the following line starts with a digit or a minus sign
            elif First_Char.isdigit() or First_Char == '-':
    #The temperature range
                 vals_list = line[0:2]
                 #try:
                 str_vals = F.readline().rstrip()+F.readline().rstrip()
                 #except UnicodeDecodeError as p:
                 #   print(p)
                 #   print(species_name)
                 str_vals = str_vals.replace('D','e')
                 str_vals = str_vals.replace('d','e')
    #Variable for the column spacing of the file
                 step = 16
    #List of coefficients
                 coef_list = [str_vals[i:i+step] for i in \
                              range(0,len(str_vals),step)]
    #Final list with Temp range and Coefficients
                 vals_list.extend(coef_list)
    #Try and except statements are for debugging. Ie even if the error occurs
    #   the program will still run instead of crashing
                 #try:
    #Reshaping the array into column array
                 vals = np.array(vals_list,dtype=float)\
                 .reshape(len(vals_list),1)
                 #except ValueError as e:
                     #print(e)
                     #print(vals_list)
                 if species_name in NEWNASA_Species:
    #Stacking the data for each temperature range into that species
                     NEWNASA_Species[species_name] = \
                     np.hstack([NEWNASA_Species[species_name],vals])
                 else:
    #Otherwise just save the data to the species
                     NEWNASA_Species[species_name] = vals
    return NEWNASA_Species, NEWNASA_SPECIES_NAMES
