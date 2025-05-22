#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 15:37:16 2020

@author: avbritt
"""

"""
This python module plots the initial and final states as bar charts.
"""
import numpy as np
from matplotlib import pyplot as plt

#import Statements (our modules)
import load_vectors

def Plot_outputs(n_true, final_n):
  global l_init
  #global input_file, names
  v_init, l_init, n_init, om, sms, scale_factor, names = \
    load_vectors.load_input()

  outs = []


  #Takes a given matrix and removes values from that matrix that are not
  #   assoiciated with the specified phase
  def onlyPhase(phaseNumber, old):
    new = []
    #For each species, if the species phase is the correct phase, add it to
    #   the new matrix (thus the new matrix will only contain species of
    #   the correct phase)
    for temp in range(len(old)):
        if l_init[temp] == phaseNumber:
    #The new version of the matrix with only the given phase
            new = np.append(new,old[temp])
    return new


  #Create figure with all the species on the same axis
  fig, ax = plt.subplots(1,1,figsize=[20,7])

  index = np.arange(len(n_true))
  bar_width = 0.35

  opacity = 0.5
  error_config = {'ecolor': '0.3'}
  ax.set_yscale('log')

  rects1 = ax.bar(index, n_true, bar_width,
                alpha=opacity, color='b',label='observed',bottom=10**-10)

  rects2 = ax.bar(index + bar_width, final_n[:len(n_true)], bar_width,bottom=10**-10,
                alpha=opacity, color='r',
                label='equilibrium')

  ax.set_ylabel('Mixing Ratio')
  #ax.set_ylabel('Scores')
  ax.set_xticks(index + bar_width / 2)
  ax.set_xticklabels(names)
  ax.legend()
  ax.set_ylim(10**-10,1)
  ax.grid()
  fig.tight_layout()
  plt.show()


  return outs
