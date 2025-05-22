#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 17:51:14 2020

@author: avbritt
"""

"""
This module takes a cell array of elements and the temperature in Kelvin and
returns those elements in "standard form" for calculating deltaG
"""
import numpy as np
import pandas as pd

import Parsing as Pars

def convert2standard(data_array,T):
    #For each element in the cell array of elements, find the letter that
    #they match.
    #print(data_array)

    #AVY each element is a dictionary entry with the element name as the key
    #   and the number of atoms as the list of data
    #Converting the dictionary into a pandas dataframe
    data_array = pd.DataFrame.from_dict(data_array)

    #If the element is diatomic, add a two afterward.
    diatomic_names = ['Br', 'I', 'N', 'Cl', 'H', 'O', 'F']
    for names in diatomic_names:
        # if names in data_array:
        if names in data_array:
            data_array[names] = [2]
            data_array = data_array.rename(columns={names:names+str('2')})
    #If the element is carbon, turn it into graphite, the standard form of
    #   carbon
    #   So, we convert the old column name from C to C(gr)
    if 'C' in data_array:
        # data_array['C(gr)'] =
        data_array = data_array.rename(columns={'C':'C(gr)'})
            #print(data_array)
        #If the element is sulfur, turn it into the reference state for sulfur
        #   at that temp.
    if 'S' in data_array:
        if T < 368.3007:
            # data_array['S(a)'] = data_array.pop('S')]
            data_array = data_array.rename(columns={'S':'S(a)'})
        elif T < 388.3607:
            # data_array['S(b)'] = data_array.pop('S')
            data_array = data_array.rename(columns={'S':'S(b)'})

        else:
            # data_array['S(L)'] = data_array.pop('S')
            data_array = data_array.rename(columns={'S':'S(L)'})

    #print(data_array)

    return data_array
