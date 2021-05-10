#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 21:05:13 2021

@author: eleonora
"""

import settings
#import force
#import velocity_verlet2
import matplotlib.pyplot as plt
import numpy as np
import printing

max_r = int(settings.t_p/20)


def energy_a():
    K_inst = np.loadtxt("./es_h/kinetic_energy_production"+ printing.names + ".txt", max_rows=max_r)
    U_inst = np.loadtxt("./es_h/potential_energy_production"+ printing.names + ".txt", max_rows=max_r)
    Etot_inst = np.loadtxt("./es_h/total_energy_production"+ printing.names + ".txt", max_rows=max_r)
    temp = np.loadtxt("./es_h/temperature_eq"+ printing.names + ".txt", max_rows=max_r)
    
    mean_K = np.mean(K_inst)
    var_K = np.var(K_inst)
    stdev_K = np.std(K_inst)
    
    mean_U = np.mean(U_inst)
    var_U = np.var(U_inst)
    stdev_U = np.std(U_inst)
    
    mean_E = np.mean(Etot_inst)
    var_E = np.var(Etot_inst)
    stdev_E = np.std(Etot_inst)
    
    ##plot energy
    x = np.arange(max_r)
    plt.figure(figsize=(12,6))
    plt.plot(x, K_inst)
    plt.hlines(mean_K, 0, len(x),color="red", label="mean")
    plt.xlabel('Number of time steps')
    plt.ylabel('K [reduced units]')
    plt.title("Kinetic energy as a function of the numeber of steps")
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(12,6))
    plt.plot(x, U_inst)
    plt.hlines(mean_U,0, len(x), color="red", label = "mean")
    plt.xlabel('Number of time steps')
    plt.ylabel('U [reduced units]')
    plt.title("Potential energy as a function of the number of steps")
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(12,6))
    plt.plot(x, Etot_inst)
    plt.xlabel('Number of time steps')
    plt.ylabel('Total energy U+K [reduced units]')
    plt.title("Total energy as a function of the number of steps")
    #plt.ylim(153.813, 153.814)
    plt.hlines(mean_E, 0, len(x), color="red",label="mean")
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(12,6))
    plt.plot(x, temp)
    plt.xlabel('Number of time steps')
    plt.ylabel('Temperature [reduced units]')
    plt.title("Temperature as a function of the number of steps")
    #plt.ylim(153.813, 153.814)
    plt.hlines(2.0, 0, len(x), color="red",label="desired temperature")
    plt.legend()
    plt.show()
    
    
    return mean_K, var_K, stdev_K, mean_U, var_U, stdev_U, mean_E, var_E, stdev_E
    
def vel_a():
    vx2_inst = np.loadtxt("./es_h/vx2_production"+ printing.names + ".txt", max_rows=max_r)
    vy2_inst = np.loadtxt("./es_h/vy2_production"+ printing.names + ".txt", max_rows=max_r)
    vz2_inst = np.loadtxt("./es_h/vz2_production"+ printing.names + ".txt", max_rows=max_r)
    
    mean_vx2 = np.mean(vx2_inst)
    var_vx2 = np.var(vx2_inst)
    stdev_vx2 = np.std(vx2_inst)
    
    mean_vy2 = np.mean(vy2_inst)
    var_vy2 = np.var(vy2_inst)
    stdev_vy2 = np.std(vy2_inst)
    
    mean_vz2 = np.mean(vz2_inst)
    var_vz2 = np.var(vz2_inst)
    stdev_vz2 = np.std(vz2_inst)
    
    
    return mean_vx2, var_vx2, stdev_vx2, mean_vy2, var_vy2, stdev_vy2, mean_vz2, var_vz2, stdev_vz2    
    