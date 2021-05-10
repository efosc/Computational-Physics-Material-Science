#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:23:43 2021

@author: eleonora
"""

import settings
import force
import velocity_verlet2
import matplotlib.pyplot as plt
import printing
import analysis
import numpy as np 

import time
start_time = time.time()


#arrays with the initial positions
x_old,y_old,z_old = settings.initial_pos()

#arrays with the initial velocities
vx_old, vy_old, vz_old = settings.initial_velocities()

#arrays with the initial forces
fx_old,fy_old,fz_old, U = force.get_forces2(x_old,y_old,z_old, 0)


#equilibration run
itime = 0
while itime < settings.t_eq:
    print(itime)
    #rescale velocities every 10 steps:
    if(itime%10==0):
        print("rescale")
        vx_old,vy_old,vz_old , K = settings.rescale_vel(vx_old,vy_old,vz_old)
    #move
    x_new, y_new, z_new, vx_new, vy_new, vz_new, fx_new, fy_new, fz_new, U = velocity_verlet2.move_vel_verl(x_old, y_old, z_old, vx_old, vy_old, vz_old,fx_old, fy_old, fz_old, itime)
    #printing
    if(itime%10==0):
        print(K,U,U+K)
        printing.write_eq(K,U)
        
    ##get ready for the next step
    x_old,y_old,z_old = x_new,y_new,z_new
    vx_old,vy_old,vz_old = vx_new,vy_new,vz_new
    fx_old,fy_old,fz_old = fx_new,fy_new,fz_new
    itime+=1 
    
#plot positions
fig1 = plt.figure(figsize=(12,9))
ax1 = fig1.add_subplot(111, projection='3d')
for i in range(settings.N):
    ax1.scatter(x_old[i], y_old[i], z_old[i], c="C0", s=20)
ax1.set_xlabel(f'x [$\sigma$ units]')
ax1.set_ylabel(f'y [$\sigma$ units]')
ax1.set_zlabel(f'z [$\sigma$ units]')
ax1.set_title("Particles' positions after equilibration")
plt.show()

#production run
itime = 0
while itime < settings.t_p:
    print(itime)
    #move
    x_new, y_new, z_new, vx_new, vy_new, vz_new, fx_new, fy_new, fz_new, U = velocity_verlet2.move_vel_verl(x_old, y_old, z_old, vx_old, vy_old, vz_old,fx_old, fy_old, fz_old, itime)
#printing
    if(itime%20==0):
        #print(K,U,U+K)
        vx2 = vx_new.dot(vx_new)
        vy2 = vy_new.dot(vy_new)
        vz2 = vz_new.dot(vz_new)
        K = 0.5*settings.m*(vx2 + vy2 + vz2)
        #printing.write_prod(vx2,vy2,vz2,K,U)
        
    ##get ready for the next step
    x_old,y_old,z_old = x_new,y_new,z_new
    vx_old,vy_old,vz_old = vx_new,vy_new,vz_new
    fx_old,fy_old,fz_old = fx_new,fy_new,fz_new
    itime+=1    

printing.close_files()


#plot positions
fig1 = plt.figure(figsize=(12,9))
ax1 = fig1.add_subplot(111, projection='3d')
for i in range(settings.N):
    ax1.scatter(x_old[i], y_old[i], z_old[i], c="C0", s=20)
ax1.set_xlabel(f'x [$\sigma$ units]')
ax1.set_ylabel(f'y [$\sigma$ units]')
ax1.set_zlabel(f'z [$\sigma$ units]')
ax1.set_title("Particles' final positions")
plt.show()

###calculate mean and standard deviation

printing.print_ave()

## point h
vel_values = vx_old**2 + vy_old**2 + vz_old**2
vel = open("final_vel_10.txt", "w")
vel_ok = open("final_vel2_10.txt", "w")
for i in range(len(vx_old)):
    vel.write(str(vx_old[i]) +"  " + str(vy_old[i]) + "  " + str(vz_old[i]) + "\n")
    vel_ok.write(str(vel_values[i]) + "\n")
  

plt.hist(vel_values, bins=np.linspace(0., max(vel_values), num=20), density=True)
plt.xlabel(f'v^2[reduced units]')
plt.ylabel("frequency")
plt.title("Velocities ditribution")

print("--- %s seconds ---" % (time.time() - start_time))

    

