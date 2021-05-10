#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 14:30:01 2021

@author: eleonora
"""

import settings
import force
import settings
import numpy as np

def move_vel_verl (x_old, y_old, z_old, vx_old, vy_old, vz_old,fx_old, fy_old, fz_old, itime):
    ## new positions: r(t+dt) = r(t) + v(t)*dt + dt^^2*f(t)/(2m)
    x_new = np.mod(x_old + settings.dt*vx_old + (settings.dt*settings.dt*fx_old)/(2.*settings.m), settings.L)
    y_new = np.mod(y_old + settings.dt*vy_old + (settings.dt*settings.dt*fy_old)/(2.*settings.m),settings.L)
    z_new = np.mod(z_old + settings.dt*vz_old + (settings.dt*settings.dt*fz_old)/(2.*settings.m),settings.L)
    # evalaute new forces
    fx_new, fy_new, fz_new,U = force.get_forces2(x_new,y_new,z_new, itime)
    ## new velocities: v(t+dt) = v(t) + dt*(f(t) + f(t+dt))/(2m)
    vx_new = vx_old + settings.dt*(fx_old + fx_new)/(2.*settings.m)
    vy_new = vy_old + settings.dt*(fy_old + fy_new)/(2.*settings.m)
    vz_new = vz_old + settings.dt*(fz_old + fz_new)/(2.*settings.m)
    return x_new, y_new, z_new, vx_new, vy_new, vz_new, fx_new, fy_new, fz_new, U
    
    
