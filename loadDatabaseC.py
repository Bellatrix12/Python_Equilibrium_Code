#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:45:06 2019

@author: avbritt
"""

"""
This program loads the aqueous phase database in the file 'F' into a dictionary
'aqueous_species'
"""

#Grabs whatever directory you are in (make sure you are in the top
#level directory of Multiphase_Equilibrium_Code)
import os
dirpathdefault = os.getcwd()


def loadDatabaseC():
    #The name of the file holding the data, which is the same as the original
    #   except for the fact that CO is renamed to find it more easily
    with open('sprons96_edited2.dat','r+') as F:
    #We need to store the species names, Temps, and coefficients associated
    #   with all the aqueous species in this data file.
    #A while loop for readline() to skip the intro comments
    #   until we get to the aqueous species section
        while F.readline().strip() != 'aqueous species':
            continue
    #Skips the line full of *'s in the file
        _=F.readline()
    #creates an empty dictionary of species
        aqueous_species = {}
    #reads the line and strips it of white space on both sides
        line = F.readline().strip()
    #While the line is not empty
        while line:
    #.split turns the line into a list of items
    #if the list has two items store the second item as the species
    #   name (this will be our dictionary key)
            if len(line.split()) == 2:
                species = line.split()[1]
    #if the line list only has one item, store that item in species
            elif len(line.split()) == 1:
                species = line
    #Bit of code to skip the next two lines
            _=F.readline()
            _=F.readline()
    #Read in the lines that have the Temperature and coefficient values,
    #   strip the lines of whitespace, and turn the line into a list
    #   of values. Once we have read in all the values for a species we use
    #   the extend() function to save all the values into one list
            vals = F.readline().strip().split()
            vals.extend(F.readline().strip().split())
            vals.extend(F.readline().strip().split())
    #Save that list of values into the dictionary key associated with the
    #   species
    #Some of the species are duplicated in the database but have different
    #   values for the Temp. and coeff. values. So the if statement below
    #   checks for duplicates and makes sure no information is overwritten
            if species in aqueous_species:
                aqueous_species[species].append(vals)
                #print(species)
            else:
                aqueous_species[species] = vals
    #Strip the next line of whitespace so that we may prepare for the next
    #   line in the while loop
            line = F.readline().strip()
    #Return the database
    return aqueous_species
