#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:09:29 2020

@author: avbritt
"""

"""
This module takes in a cell array of species names and creates a coefficient
matrix for it.
NOTE: the top row has all the elements in it.
Precondition: names must have at least one species.
"""

import numpy as np
import pandas as pd
import re

def makeAs(names):

    #~*~*~*~**~*~*~*~**~*~FUNCTIONS SECTION ~*~*~*~**~*~*~*~**~*~*~**~*~*~*~**
    #This function is to parse species names like CH4 into is constituent
    #   elements (C 1 H 4) and save it to a dictionary with the name of
    #   the element as the dictionary key.
    def parse(species):
    #initializing an empty dictionary
        d = {}

    #Expression to strip off any any characters surrounded by '()'
    #   to take care of species like, for example, charged species
    #   like CO3(-2), aqueous species like NH3(0), etc.
        species = re.sub('\(.+\)','',species)
    #Strip off any non-alphanumeric characters in the species
        #species = re.sub(r'\W+','',species)
    #Strip off _L from H2O_L species
        if species == 'H2O_L':
            species = re.sub('[_L]','',species)
    #Initializing empty strings for the element name and the number (num)
    #   of atoms for that element
        name = ''
        num = ''
    #A for loop to loop over each character in the species name
        for char in species:
    #Check to see if the character is uppercase. If it is, clear name and num
    #   variable and save that character to name variable. On the next
    #   iteration then, if name and num are filled we can save them to the
    #   dictionary and if num is not filled on this iteration, then the
    #   number of atoms for that element is assumed to be 1 and we save
    #   name and num to the dictionary.
            if char.isupper():
                    if name:
                        if not num:
                            num = '1'
                        d[name] = [int(num)]
                    name = ''
                    num = ''
                    name+=char
    #Check to see if the character is lower case. If it is, save that
    #  character to the name variable
            elif char.islower():
                    #print('char was lower')
                name+=char
    #Check to see if the character is numeric. If so, save that number
    #  to the num variable
            elif char.isdigit():
                num+=char
    #On the first iteration, is num is empty then we assume the number of
    #   atoms to be 1 and we save name and num to the dictionary
        if not num:
            num = '1'
        d[name] = [int(num)]
   #Clear name and num variables in preparation for the next iteration
        name = ''
        num = ''
   #     print(d)
        return d



    #~*~*~**~*~*~*~*~**~*~*~**~*~*~*~**~*~*~*~***~*~*~**~*~*~**~*~*~*~**~*~*

    #Loop over each species and call the parse function to parse each species
    #   into its constituent elements
    full_dict = {}
   # for species in names:
    #Call to the parse function with takes in one species at a time and
    #   separates it into its constituent elements and returns a
    #   single species dictionary
    d = parse(names)
    #print(species,d)
    #for loop for the elements within the species dictionary
    for dict_species in d:
        #print(np.size(d[dict_species]))
        #print(dict_species)
        #print(d[dict_species])
        #full_dict[dict_species] = d[dict_species]
    #full_dict becomes a dictionary of a singular species divided into the
    #   elements it's composed of
        if dict_species in full_dict:
            full_dict[dict_species].append(d[dict_species])
          #full_dict.update(d)
        else:
            full_dict[dict_species]= d[dict_species]

    return full_dict
