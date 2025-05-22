#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module loads the input text file containing species names,
abundances, charges, and species type, and stores them as arrays which are then
passed back to the main code
"""
#Grabs whatever directory you are in (make sure you are in the top
#level directory of Multiphase_Equilibrium_Code)
import os
dirpathdefault = os.getcwd()

import numpy as np
import pandas



def load_input():
  #Choose input file (uncomment the input file you wish to use, or create your own)
  #Defaults are gas phase templates for Modern Earth, and Mars

  #Template files created by AVY
  #input_file = 'inputs_Mars.txt'
  input_file = 'inputs_ModernEarth_gas_only.txt'

  #number of species included in the template (ex defaults:
  # Modern Earth gas only system (15), Mars (14),
  n_spec = 15

  #Load the data from the text input file that was chosen
  All_inputs = pandas.read_csv(input_file, names=list(range(n_spec)))

  #First line is the charge vector (v_input-->v_init)
  v_input = All_inputs[0][0]
  #splits the variable on whitespace so we have a list of multiple items
  v_init = v_input.split()
  #Changes the datatype of each list item to floats
  v_init = [int(i) for i in v_init]

  #Second line is the phase vector (l_input-->l_init)
  l_input = All_inputs[0][1]
  #splits the variable on whitespace so we have a list of multiple items
  l_init = l_input.split()
  #Changes the datatype of each list item to floats
  l_init = [int(i) for i in l_init]

  #Fourth line is the abundance vector (n_input-->n_init)
  n_input = All_inputs[0][3]
  #splits the variable on whitespace so we have a list of multiple items
  n_init = n_input.split()
  #Changes the datatype of each list item to floats
  n_init = [float(i) for i in n_init]

  #Fifth line contains various scalars (scalar_input-->scalars_in)
  # np.asarray evaluates the pandas dataframe as an array so I can
  # grab the 5th row
  scalar_input = np.asarray(All_inputs)[4,:]
  #Making a an array of only the first three elements since the rest are empty
  scalars_in = scalar_input[0:3]
  #Changes the datatype of each list item to floats
  scalars_in = [float(i) for i in scalars_in]

  #Third line contains the names of the species (name_input)
  #   (1x200 character list for Archean_max)
  name_input = np.asarray(All_inputs)[2,:]
  name_input = name_input.tolist()


  #scalar_in is an array (list) of the number of ocean masses (om),
  #   salinity scalar (sms) (which is no longer used), and scale factor
  #   (scale_factor)
  om = scalars_in[0]
  sms = scalars_in[1]
  scale_factor = scalars_in[2]

  #Multiply liquid water abundance by number of ocean masses
  #   n_init(1)=om*n_init(1)
  n_init[0]=om*n_init[0]

  return v_init, l_init, n_init, om, sms, scale_factor, name_input
