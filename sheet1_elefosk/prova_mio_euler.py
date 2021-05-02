#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 10:54:54 2021

@author: eleonora
"""
import numpy as np
import matplotlib.pyplot as plt

#global variables
n_planets = 9
sun_mass = np.loadtxt("mass.dat", skiprows=3, usecols=0, max_rows=1)*1E29
masses = np.loadtxt("mass.dat", skiprows=4, usecols=0)*1E29

AU = 1.49597870691E+11
AUday = 1.49597870691E+11/86400.
G = 6.67408E-11/(AUday**2)/AU

def get_forces(x,y,z): ## fij = G*mi*mj*(xj-xi)/(r^3)
    fx = np.zeros(9,dtype=float)
    fy = np.zeros(9,dtype=float)
    fz = np.zeros(9,dtype=float)
    i=0
    while i< n_planets:
        ## consider the ith planet 
        xi = x[i]
        yi = y[i]
        zi = z[i]
        mi = masses[i]
        ##sun
        fxi = -(G*mi*sun_mass*(xi))/((np.sqrt(xi**2 + yi**2 +zi**2))**3)
        fyi = -(G*mi*sun_mass*(yi))/((np.sqrt(xi**2 + yi**2 +zi**2))**3)
        fzi = -(G*mi*sun_mass*(zi))/((np.sqrt(xi**2 + yi**2 +zi**2))**3)
        ## all other planets:
        j=0
        while j<n_planets:
            if(i!=j):
                xj = x[j]
                yj = y[j]
                zj = z[j]
                mj = masses[j]
                fxi += (G*mi*mj*(xj-xi))/((np.sqrt((xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2))**3)
                fyi += (G*mi*mj*(yj-yi))/((np.sqrt((xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2))**3)
                fzi += (G*mi*mj*(zj-zi))/((np.sqrt((xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2))**3)
            j+=1 
        fx[i], fy[i],fz[i] = fxi,fyi,fzi
        i+=1
    
    return fx,fy,fz

##potential energy of the Solar System 
def get_U(x,y,z):
    i=0
    U=0
    while i< n_planets:
        ## consider the ith planet 
        xi = x[i]
        yi = y[i]
        zi = z[i]
        mi = masses[i]
        ##sun
        U += -(G*mi*sun_mass)/(np.sqrt(xi**2 + yi**2 +zi**2))
        ## all other planets:
        j=0
        while j<n_planets:
            if(i!=j):
                xj = x[j]
                yj = y[j]
                zj = z[j]
                mj = masses[j]
                U += -(G*mi*mj)/(np.sqrt((xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2))
            j+=1 
        i+=1
    
    return U

## Kinetic energy of the Solar System
def get_K(vx,vy,vz):
    i=0
    K=0
    while i< n_planets:
        ## consider the ith planet 
        vxi = vx[i]
        vyi = vy[i]
        vzi = vz[i]
        mi = masses[i]
        K += (mi*(vxi)**2 + (vyi)**2 +(vzi)**2)/2.
        i+=1
    
    return K


#simulation conditions
imax = 73000-1     # number of iterations
dt = 0.5    # time step
itime=0

## prepare array for total energies throughout the simulation
total_energy = np.zeros(int((imax+1)/100))

## save Pluto's position to plot trajectory
pluto_x = np.zeros(int((imax+1)/100))
pluto_y = np.zeros(int((imax+1)/100))
pluto_z = np.zeros(int((imax+1)/100))

##arrays with the "old" position each time -> at the beginning initial positions
x_old = np.loadtxt("planets.dat", skiprows=7, usecols=0)
y_old = np.loadtxt("planets.dat", skiprows=7, usecols=1)
z_old = np.loadtxt("planets.dat", skiprows=7, usecols=2)

##arrays with "old velocities for the first step- > at the beginning initial velocities
vx_old = np.loadtxt("planets.dat", skiprows=7, usecols=3)
vy_old = np.loadtxt("planets.dat", skiprows=7, usecols=4)
vz_old = np.loadtxt("planets.dat", skiprows=7, usecols=5)

##save pluto's first position
pluto_x[0], pluto_y[0], pluto_z[0] = x_old[8], y_old[8],z_old[8]
## initial energy
total_energy[0] = get_U(x_old,y_old,z_old) + get_K(vx_old,vy_old,vz_old)

while itime<imax:
    print(itime)
    ##evaluate forces from "old" positions
    fx,fy,fz = get_forces(x_old, y_old,z_old)
    ##new positions: r(t+dt) = r(t) + v(t) + f(t)*dt^^2/(2*m)
    x_new = x_old + dt*vx_old + (dt*dt*fx)/(2.*masses)
    y_new = y_old + dt*vy_old + (dt*dt*fy)/(2.*masses)
    z_new = z_old + dt*vz_old + (dt*dt*fz)/(2.*masses)
    ##new velocities: v(t+dt) = v(t) + dt*f/m
    vx_new = vx_old + (dt*fx)/(masses)
    vy_new = vy_old + (dt*fy)/(masses)
    vz_new = vz_old + (dt*fz)/(masses)
    ##save energy and plutos positions (1 every 100)
    if((itime+1)%100 == 0):
        pluto_x[int((itime+1)/100)],pluto_y[int((itime+1)/100)],pluto_z[int((itime+1)/100)] = x_new[8],y_new[8],z_new[8]
        total_energy[int((itime+1)/100)] = get_U(x_new,y_new,z_new) + get_K(vx_new,vy_new,vz_new)
    ##get ready for the new step
    x_old,y_old,z_old = x_new,y_new,z_new
    vx_old,vy_old,vz_old = vx_new,vy_new,vz_new 
    itime+=1
    
#plot plutos trajectory
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(pluto_x,pluto_y,pluto_z)
ax.set_xlabel('x [a.u.]')
ax.set_ylabel('y [a.u.]')
ax.set_zlabel('z [a.u.]')
ax.scatter(0,0,0,c="red", s=50, label="Sun") #sun position
ax.legend()
ax.set_title("Pluto's trajectory in 300 years, Euler algorithm")
plt.show()

##plot energy
x = np.arange(len(total_energy))*100*43200  ## half day (1 measurement=100 steps) -> seconds
plt.figure(figsize=(12,6))
plt.plot(x, total_energy*AUday**2)
plt.xlabel('Time [s]')
plt.ylabel('Total energy U+K [J]')
plt.title("Total energy 100 years, Euler's algorithm")
plt.show()