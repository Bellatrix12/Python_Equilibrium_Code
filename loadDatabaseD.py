#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:59:02 2019

@author: avbritt
"""

"""
This program loads the 4th database in the file 'F' into a pandas dataframe'
"""
#Grabs whatever directory you are in (make sure you are in the top
#level directory of Multiphase_Equilibrium_Code)
import os
dirpathdefault = os.getcwd()

import numpy as np
import pandas

def loadDatabaseD():
    #The name of the file holding the data. This file has the
    #   undefined character in the first line removed in order to read the file
    F = pandas.read_csv('AQ_coefficients1.txt', header=None,delimiter='\t', \
                        encoding='UTF-16BE')
    #F = np.genfromtxt('AQ_coefficients1.txt',dtype=str)

    """
    We set the column names as
    [species_name, a, b] and then set the row index to the species names
    That way we can access the a and b value of each species using the
    Database.loc['species_name']
    """
    F.columns = ['species_name','a','b']
    F = F.set_index('species_name')

    return F
