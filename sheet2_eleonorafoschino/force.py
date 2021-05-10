#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 19:25:01 2021

@author: eleonora
"""
import numpy as np
import settings 
## x,y,z are vectors with the positions of all the particles
r_cut = 2.5*settings.sigma
U_rcut = 4*settings.epsilon*(settings.sigma**12/(r_cut**12) - settings.sigma**6/(r_cut**6))

    
def pbc(x):
    return x - np.rint(x/settings.L)*settings.L
   
def get_forces(x,y,z): ## calculate force between two particles only if r<r_cut
    uno = np.ones(shape=(settings.N,settings.N))
    X = np.transpose(x*uno)
    Y = np.transpose(y*uno)
    Z = np.transpose(z*uno)
    
    r_x = pbc(X - np.transpose(X)) 
    r_y = pbc(Y - np.transpose(Y)) 
    r_z = pbc(Z - np.transpose(Z))
    r = np.sqrt(r_x*r_x + r_y*r_y + r_z*r_z) - r_cut #shift potential
    r_rec = np.reciprocal(r) 
    r_rec[r_rec > 0] = 0. ## set to 0  if distance > r_cut
    f_allx = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_x
    fx = f_allx.sum(axis=0)
    f_ally = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_y
    fy = f_ally.sum(axis=0)
    f_allz = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_z
    fz = f_allz.sum(axis=0)
    return fx, fy, fz


def get_forces2(x,y,z, itime):  ##shifts potential energy and calculate U
    
    uno = np.ones(shape=(settings.N,settings.N))
    X = np.transpose(x*uno)
    Y = np.transpose(y*uno)
    Z = np.transpose(z*uno)
    
    r_x = pbc(X - np.transpose(X)) 
    r_y = pbc(Y - np.transpose(Y)) 
    r_z = pbc(Z - np.transpose(Z))
    r = np.sqrt(r_x*r_x + r_y*r_y + r_z*r_z) +np.eye(settings.N) 
    r_rec = np.reciprocal(r) -np.eye(settings.N)
    f_allx = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_x
    fx = f_allx.sum(axis=0)
    f_ally = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_y
    fy = f_ally.sum(axis=0)
    f_allz = 4.*settings.epsilon*(-12.*settings.sigma**12*r_rec**14 + 6.*settings.sigma**6*r_rec**8)*r_z
    fz = f_allz.sum(axis=0)
    # calculate U every 10 steps
    if(itime%10 == 0):
        U_all = 4.*settings.epsilon*(settings.sigma**12*r_rec**12 - settings.sigma**6*r_rec**6) + U_rcut
        U_tot = U_all.sum()/2.
    else:
        U_tot = 0.
    return fx, fy, fz, U_tot

