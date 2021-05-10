#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  9 15:23:10 2021

@author: eleonora
"""

import settings
import force
import velocity_verlet2
import matplotlib.pyplot as plt
import analysis


#open files for writing
names = "NO"

ek_eq = open("es_h/kinetic_energy_eq" + names + ".txt", "w")
u_eq = open("es_h/potential_energy_eq"+ names + ".txt", "w")
etot_eq = open("es_h/total_energy_eq"+ names + ".txt", "w")
    #vel_eq = open("velocities_eq.txt", "w")
    #pos_eq = open("positions_eq.txt", "w")
temp_eq = open("es_h/temperature_eq"+ names + ".txt", "w")

ek_production = open("es_h/kinetic_energy_production"+ names + ".txt", "w")
u_production = open("es_h/potential_energy_production"+ names + ".txt", "w")
etot_production = open("es_h/total_energy_production"+ names + ".txt", "w")
#vel_eq = open("velocities_eq.txt", "w")
    #pos_eq = open("positions_eq.txt", "w")
temp_production = open("es_h/temperature_production"+ names + ".txt", "w")
vx2_production = open("es_h/vx2_production"+ names + ".txt", "w")
vy2_production = open("es_h/vy2_production"+ names + ".txt", "w")
vz2_production = open("es_h/vz2_production"+ names + ".txt", "w")
averages = open("es_h/averages" + names+".txt", "w")


def write_eq(K,U):
    ek_eq.write(str(K) + "\n")
    u_eq.write(str(U) + "\n")
    etot_eq.write(str(K + U) + "\n")
    temp_eq.write (str(2./3.*settings.kB*K/float(settings.N)) + "\n")
    

def write_prod(vx2,vy2,vz2,K,U):
    vx2_production.write(str(vx2) + "\n")
    vy2_production.write(str(vy2) + "\n")
    vz2_production.write(str(vz2) + "\n")
    ek_production.write(str(K) + "\n")
    u_production.write(str(U) + "\n")
    etot_production.write(str(K + U) + "\n")
    temp_production.write (str(2./3.*settings.kB*K/float(settings.N)) + "\n")    
    
def close_files():
    ek_eq.close()
    u_eq.close()
    etot_eq.close()
    temp_eq.close()
    ek_production.close()
    u_production.close()
    etot_production.close()
    temp_production.close()
    vx2_production.close()
    vy2_production.close()
    vz2_production.close()
    
def print_ave():

    mean_K, var_K, stdev_K, mean_U, var_U, stdev_U, mean_E, var_E, stdev_E = analysis.energy_a()
    mean_vx2, var_vx2, stdev_vx2, mean_vy2, var_vy2, stdev_vy2, mean_vz2, var_vz2, stdev_vz2 = analysis.vel_a()
    
    print("mean K, var_K, stdev_K \n")
    print(str(mean_K) + " "+ str(var_K)  + " "+ str(stdev_K) + "\n")
    print("mean U, var_U, stdev_U \n")
    print(str(mean_U) + " "+  str(var_U) + " "+ str(stdev_U) + "\n")
    print("mean E, var_E, stdev_E \n")
    print( str(mean_E) + " "+ str(var_E) + " "+ str(stdev_E) + "\n")

    print("mean_vx2, var_vx2, stdev_vx2 \n")
    print(str(mean_vx2) + " "+  str(var_vx2) + " "+ str(stdev_vx2) + "\n")
    print("mean_vy2, var_vy2, stdev_vy2 \n")
    print(str(mean_vy2) + " "+ str(var_vy2) + " "+ str(stdev_vy2) + "\n")
    print("mean_vz2, var_vz2, stdev_vz2 \n")
    print(str(mean_vz2) + " "+ str(var_vz2) + " "+ str(stdev_vz2) + "\n")
    
    averages.write("mean K, var_K, stdev_K \n\n")
    averages.write(str(mean_K) + " "+ str(var_K)  + " "+ str(stdev_K) + "\n")
    averages.write("mean U, var_U, stdev_U \n")
    averages.write(str(mean_U) + " "+  str(var_U) + " "+ str(stdev_U) + "\n")
    averages.write("mean E, var_E, stdev_E \n")
    averages.write( str(mean_E) + " "+ str(var_E) + " "+ str(stdev_E) + "\n")

    averages.write("mean_vx2, var_vx2, stdev_vx2 \n")
    averages.write(str(mean_vx2) + " "+  str(var_vx2) + " "+ str(stdev_vx2) + "\n")
    averages.write("mean_vy2, var_vy2, stdev_vy2 \n")
    averages.write(str(mean_vy2) + " "+ str(var_vy2) + " "+ str(stdev_vy2) + "\n")
    averages.write("mean_vz2, var_vz2, stdev_vz2 \n")
    averages.write(str(mean_vz2) + " "+ str(var_vz2) + " "+ str(stdev_vz2) + "\n")
