#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 20:37:29 2020

@author: avbritt
"""

"""
Module to parse a species into its constituent elements.
The return output is a variable with the list of elements within that species
"""
import numpy as np
import re

def parser(species):
 #initializing an empty list
     d = []

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
 #   variable and save that character to name variable.
         if char.isupper():
                 if name:
                     if not num:
                         num = '1'
                     d.append(name)
#                     print(d)
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
 #On the first iteration, if num is empty then we assume the number of
 #   atoms to be 1
     if not num:
         num = '1'
     d.append(name)
#     print(d)
#Clear name and num variables in preparation for the next iteration
     name = ''
     num = ''
#     print(d)
     return d
