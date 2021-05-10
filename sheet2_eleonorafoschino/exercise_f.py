#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  9 20:28:07 2021

@author: eleonora
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

#prova fit
x = np.array([16, 36,64,100 ])
times = np.array([11.156896114349365, 82.19901585578918, 477.8447768688202, 1851.400703907013])


def tofit_f (N ,alpha):
    return N**alpha

params = scipy.optimize.curve_fit(tofit_f, xdata=x, ydata=times)

print(params)

#grafico
plt.figure(figsize=(12,9))
plt.scatter(x,times,label="simulation time values", s=50)
plt.plot(x, tofit_f(x,params[0]), label="Fit", c="C1")
plt.title("Simulation time as a function of N", size=20)
plt.xlabel('N',size=16)
plt.ylabel('simulation time ',size=16)
plt.legend(loc="lower right", fontsize=16)
plt.grid(True)