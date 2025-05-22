#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 16:22:30 2020

@author: avbritt
"""

"""
Module that calculates the dielectric constant for water using the data from
Owin 1961
"""
import numpy as np

def dielectric(T,P):
    t = T-273.15
    p = P-1
    a = [-22.5713, -0.32066, -0.00028568, 0.001182, 0.000027895, \
         -0.00000001476, 2300.64, -0.13476]
    D0 = 4.476150
    di = np.exp((-10**6*D0+2*a[0]*p+2*a[1]*p*t+2*a[2]*p*t**2+2*a[3]*p**2+2\
                 *a[4]*p**2*t+2*a[5]*p**2*t**2+2*a[6]*t+2*a[7]*t**2)/(-(10**6)))
    return di
