#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 18:01:22 2020

@author: avbritt
"""

"""
This module loads the database of Pitzer binary interaction parameters from
their database file (pitzer_coef.txt) and saves them in a dictionary
"""

#Grabs whatever directory you are in (make sure you are in the top
#level directory of Multiphase_Equilibrium_Code)
import os
dirpathdefault = os.getcwd()


def load_Pitzer():
    #Create an empty dictionary of Pitzer parameters
    Pitzer_params = {}
    #Read in the name of the file holding the data
    with open('pitzer_coef.txt','r+') as F:
    #reads the line and strips it of white space on both sides
        line = F.readline().strip()
    #While the line is not empty
        while line:
            if len(line.split()) == 1:
                cell_name = line.split()[0][1:]
            elif len(line.split()) >1:
                vals = line.split()
                #print(vals)
                #print(len(vals))
            #vals.extend(F.readline().strip().split())
                if cell_name in Pitzer_params:
                    Pitzer_params[cell_name].append(vals)
                else:
                    Pitzer_params[cell_name] = [vals]
    #Strip the next line of whitespace so that we may prepare fir the next
    #   line in the while loop
            line = F.readline().strip()


    return Pitzer_params
