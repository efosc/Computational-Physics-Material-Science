#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 19:19:49 2021

@author: eleonora
"""

import numpy as np
from itertools import product
import matplotlib.pyplot as plt

#global variables
kB = 1.
sigma = 1.  ## per ora esprimo tutto in unità di sigma
T = 2.
epsilon = 1.
dt=2E-4
t_eq =2000
t_p = 2000
m = 1.
imax = 100

##box 
n=10
N = n**3 ## number of particles
L = n*sigma*np.cbrt(10)
a_lat = L/n
p = L/sigma
rho = N/(L**3)


##initial positions
def initial_pos():
    ##arrays with the particles' positions
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    n_part = 0
    for i,j,k in product(list(np.arange(n)), repeat=3):
        x[n_part] = i*a_lat
        y[n_part] = j*a_lat
        z[n_part] = k*a_lat
        n_part += 1
        
    """#plot positions
    fig1 = plt.figure(figsize=(12,9))
    ax1 = fig1.add_subplot(111, projection='3d')
    for i in range(N):
        ax1.scatter(x[i], y[i], z[i], c="C0", s=20)
    ax1.set_xlabel(f'x [$\sigma$ units]')
    ax1.set_ylabel(f'y [$\sigma$ units]')
    ax1.set_zlabel(f'z [$\sigma$ units]')
    ax1.set_title("Particles' initial positions")
    plt.show()"""
    
    return(x,y,z)
    
def initial_velocities():
    #random value between -1 and 1 and then shift 
    vx0 = np.random.rand(N)*2. -1.
    vx0 = vx0 - np.mean(vx0)
    vy0 = np.random.rand(N)*2. -1.
    vy0 = vy0 - np.mean(vy0)
    vz0 = np.random.rand(N)*2. -1.
    vz0 = vz0 - np.mean(vz0)
    
    """#calculate current average velocity and rescale
    v_ave = (vx0.dot(vx0) + vy0.dot(vy0) + vz0.dot(vz0))/float(N)
    vx0 = np.sqrt(T*3/v_ave)*vx0
    vy0 = np.sqrt(T*3/v_ave)*vy0
    vz0 = np.sqrt(T*3/v_ave)*vz0"""
    
    return vx0, vy0, vz0

def rescale_vel(vx,vy,vz):
    #calculate current average velocity and rescale
    v_ave = (vx.dot(vx) + vy.dot(vy) + vz.dot(vz))/float(N)
    vx_res = np.sqrt(T*3/v_ave)*vx
    vy_res = np.sqrt(T*3/v_ave)*vy
    vz_res = np.sqrt(T*3/v_ave)*vz
    K = 0.5 * m*v_ave*float(N)    
    return vx_res,vy_res,vz_res, K
    
    
    