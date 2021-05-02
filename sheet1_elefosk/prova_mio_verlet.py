#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 09:51:18 2021

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


#simulation conditions
imax = 2190-2     # number of iterations
dt = 0.5    # time step
itime=0

mars_x = np.zeros(imax+2)
mars_y = np.zeros(imax+2)
mars_z = np.zeros(imax+2)

##arrays with the "old" position each time -> at the beginning initial positions
x_old = np.loadtxt("planets.dat", skiprows=7, usecols=0)
y_old = np.loadtxt("planets.dat", skiprows=7, usecols=1)
z_old = np.loadtxt("planets.dat", skiprows=7, usecols=2)

##arrays with iniital velocities for the first step
vx0 = np.loadtxt("planets.dat", skiprows=7, usecols=3)
vy0 = np.loadtxt("planets.dat", skiprows=7, usecols=4)
vz0 = np.loadtxt("planets.dat", skiprows=7, usecols=5)

##first simulation step
##evaluate forces
fx,fy,fz = get_forces(x_old, y_old,z_old)
## new positions: r1 = r0 + dt*v0 + dt^^2*fx/(2m)
x_mid = x_old + dt*vx0 + dt*dt*fx/(2.*masses)
y_mid = y_old + dt*vy0 + dt*dt*fy/(2.*masses)
z_mid = z_old + dt*vz0 + dt*dt*fz/(2.*masses)

##save mars' first 2 positions
mars_x[0], mars_y[0], mars_z[0] = x_old[3], y_old[3],z_old[3]
mars_x[1], mars_y[1], mars_z[1] = x_mid[3], y_mid[3],z_mid[3]

## find pluto's major semiaxis to calculate the period
max_dist2 = x_old[3]**2 + y_old[3]**2+z_old[3]**2
min_dist2 = x_old[3]**2 + y_old[3]**2+z_old[3]**2

while itime<imax:
    print(itime)
    ##evaluate forces from "middle" positions
    fx,fy,fz = get_forces(x_mid, y_mid,z_mid)
    ##new positions : r(t+dt) = 2r(t) + r(t-dt) + dt^^2*f/m
    x_new = 2*x_mid - x_old + (dt*dt*fx)/(masses)
    y_new = 2*y_mid - y_old + (dt*dt*fy)/(masses)
    z_new = 2*z_mid - z_old + (dt*dt*fz)/(masses)
    ## save mars positions
    mars_x[itime+2],mars_y[itime+2],mars_z[itime+2] = x_new[3],y_new[3],z_new[3]
    ## find mars' afelio and perielio to calculate the period
    R2 = x_new[3]**2+y_new[3]**2+z_new[3]**2
    r2 = x_new[3]**2+y_new[3]**2+z_new[3]**2
    if(R2>max_dist2): 
        max_dist2=R2
    if(r2<min_dist2): 
        min_dist2=r2
    ##get ready for the new step
    x_old,y_old,z_old = x_mid,y_mid,z_mid
    x_mid,y_mid,z_mid = x_new,y_new,z_new
    itime+=1
    
#plot Mars trajectory
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(mars_x,mars_y,mars_z)
ax.set_xlabel('x [a.u.]')
ax.set_ylabel('y [a.u.]')
ax.set_zlabel('z [a.u.]')
ax.scatter(0,0,0,c="red", s=50, label="Sun") #sun position
ax.legend()
ax.set_title("Mars' trajectory over 2 years, Verlet algorithm")
plt.show()
    
##mars orbital period according to Kepler 3rd law
major_semiaxis = (np.sqrt(R2)+np.sqrt(r2))/2.
T = np.sqrt((4.*np.pi**2/G/sun_mass)*(major_semiaxis**3))

print(T/365.)