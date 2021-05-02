#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 11:58:00 2021

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

def get_forces(x,y,z):
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
imax = 219000-1     # number of iterations
dt = 0.5    # time step
itime=0

pluto_x = np.zeros(int((imax+1)/100))
pluto_y = np.zeros(int((imax+1)/100))
pluto_z = np.zeros(int((imax+1)/100))

total_energy = np.zeros(int((imax+1)/100))

##arrays with the "old" position each time -> at the beginning initial positions
x_old = np.loadtxt("planets.dat", skiprows=7, usecols=0)
y_old = np.loadtxt("planets.dat", skiprows=7, usecols=1)
z_old = np.loadtxt("planets.dat", skiprows=7, usecols=2)

##arrays with "old" velocities at each time -> at the beginnnig initial velocities
vx_old = np.loadtxt("planets.dat", skiprows=7, usecols=3)
vy_old = np.loadtxt("planets.dat", skiprows=7, usecols=4)
vz_old = np.loadtxt("planets.dat", skiprows=7, usecols=5)

##save pluto's first position
pluto_x[0], pluto_y[0], pluto_z[0] = x_old[8], y_old[8],z_old[8]

"""pluton_x = np.zeros(imax+1)
pluton_y = np.zeros(imax+1)
pluton_z = np.zeros(imax+1)
##save pluto's first position
pluton_x[0], pluton_y[0], pluton_z[0] = x_old[8], y_old[8],z_old[8]"""

#itime=0
## find pluto's major semiaxis to calculate the period
max_dist2 = x_old[8]**2 + y_old[8]**2+z_old[8]**2
min_dist2 = x_old[8]**2 + y_old[8]**2+z_old[8]**2

while itime < imax:
    print(itime)
    ## evaluate "old" forces
    fx_old, fy_old, fz_old = get_forces(x_old,y_old,z_old)
    ## new positions
    x_new = x_old + dt*vx_old + (dt*dt*fx_old)/(2.*masses)
    y_new = y_old + dt*vy_old + (dt*dt*fy_old)/(2.*masses)
    z_new = z_old + dt*vz_old + (dt*dt*fz_old)/(2.*masses)
    # evalaute new forces
    fx_new, fy_new, fz_new = get_forces(x_new,y_new,z_new)
    ## new velocities
    vx_new = vx_old + dt*(fx_old + fx_new)/(2.*masses)
    vy_new = vy_old + dt*(fy_old + fy_new)/(2.*masses)
    vz_new = vz_old + dt*(fz_old + fz_new)/(2.*masses)
    ##save energy and plutos positions (1 every 100)
    if((itime+1)%100 == 0):
        pluto_x[int((itime+1)/100)],pluto_y[int((itime+1)/100)],pluto_z[int((itime+1)/100)] = x_new[8],y_new[8],z_new[8]
        total_energy[int((itime+1)/100)] = get_U(x_new,y_new,z_new) + get_K(vx_new,vy_new,vz_new)
    ## save the planets positions after 50 years
    ##if(itime==36500-1):
        #planets_50 = np.matrix([x_new, y_new, z_new])
    ## find pluto's afelio and perielio to calculate the period
    R2 = x_new[8]**2+y_new[8]**2+z_new[8]**2
    r2 = x_new[8]**2+y_new[8]**2+z_new[8]**2
    if(R2>max_dist2): 
        max_dist2=R2
    if(r2<min_dist2): 
        min_dist2=r2
    ##get ready for the next step
    x_old,y_old,z_old = x_new,y_new,z_new
    vx_old,vy_old,vz_old = vx_new,vy_new,vz_new
    itime+=1
 
##pluto's orbital period according to Kepler 3rd law
major_semiaxis = (np.sqrt(R2)+np.sqrt(r2))/2.
T = np.sqrt((4.*np.pi**2/G/sun_mass)*(major_semiaxis**3))

print(T/365)
    
"""#plot plutos trajectory
fig3 = plt.figure(figsize=(12,9))
ax3 = fig3.add_subplot(111, projection='3d')
ax3.plot(pluton_x,pluton_y,pluton_z)
ax3.set_xlabel('x [a.u.]')
ax3.set_ylabel('y [a.u.]')
ax3.set_zlabel('z [a.u.]')
ax3.scatter(0,0,0,c="red", s=50, label="Sun") #sun position
ax3.legend()
ax3.set_title("Pluto's trajectory in the next 100 years")
plt.show()"""
    
#planets_100 = np.matrix([x_new, y_new, z_new])

"""#plot plutos trajectory
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(pluto_x,pluto_y,pluto_z)
ax.set_xlabel('x [a.u.]')
ax.set_ylabel('y [a.u.]')
ax.set_zlabel('z [a.u.]')
ax.scatter(0,0,0,c="red", s=50, label="Sun") #sun position
ax.legend()
ax.set_title("Pluto's trajectory in 300 years")
plt.show()"""

"""#plot positions after 50 years
fig1 = plt.figure(figsize=(12,9))
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(0,0,0, c="red", s=50, label="Sun") #sun position
for i in range(9):
    ax1.scatter(planets_50[0,i], planets_50[1,i],planets_50[2,i],c="C0", s=20)
ax1.set_xlabel('x [a.u.]')
ax1.set_ylabel('y [a.u.]')
ax1.set_zlabel('z [a.u.]')
ax1.legend()
ax1.set_title("Planets' positions after 50 years")
plt.show()

#plot positions after 100 years
fig2 = plt.figure(figsize=(12,9))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(0,0,0, c="red", s=50, label="Sun") #sun position
for i in range(9):
    ax2.scatter(planets_100[0,i], planets_100[1,i],planets_100[2,i], c="C0", s=20 )
ax2.set_xlabel('x [a.u.]')
ax2.set_ylabel('y [a.u.]')
ax2.set_zlabel('z [a.u.]')
ax2.set_title("Planets' positions after 100 years")
ax2.legend()

plt.show()

##arrays with the "old" position each time -> at the beginning initial positions
x_old = np.loadtxt("planets.dat", skiprows=7, usecols=0)
y_old = np.loadtxt("planets.dat", skiprows=7, usecols=1)
z_old = np.loadtxt("planets.dat", skiprows=7, usecols=2)

##arrays with "old velocities for the first step- > at the beginning initial velocities
vx_old = np.loadtxt("planets.dat", skiprows=7, usecols=3)
vy_old = np.loadtxt("planets.dat", skiprows=7, usecols=4)
vz_old = np.loadtxt("planets.dat", skiprows=7, usecols=5)

##save pluto's first position
#pluto_x[0], pluto_y[0], pluto_z[0] = x_old[8], y_old[8],z_old[8]
## initial energy
total_energy[0] = get_U(x_old,y_old,z_old) + get_K(vx_old,vy_old,vz_old)

##plot energy
x = np.arange(len(total_energy))*100*43200
plt.figure(figsize=(12,6))
plt.plot(x, total_energy*AUday**2)
plt.xlabel('Time [s]')
plt.ylabel('Total energy U+K [J]')
plt.title("Total energy 300 years,Velocity Verlet algorithm")
plt.show()"""
