#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# This code is only meant to build the uniform distribution of points across a sphere
# This should only need to be run once per measurment sphere 
# =============================================================================

from numpy import linspace, sin, cos, arccos, pi, array

u = linspace(0,1,200)
v = linspace(1,0,200)

theta = arccos(2*v-1)
phis = 2*pi*u



def points(r,t,p):
    # Converts the spherical polars to cartesian coords
    x = r*sin(t)*cos(p)
    y = r*sin(t)*sin(p)
    z = r*cos(t)
    
    return {'x':x,'y':y,'z':z}

# Working at a radius of 1m  to ensure far field regime is achived
coords = array([[points(1e6,t_val,phi_val) for t_val in theta] for phi_val in phis])

# We need three files like this to make sure that we get different named files for the scattered fields
with open('EvalPoints.dat','w') as file:
    for cart in coords.flatten():
        file.write(f"{cart['x']} {cart['y']} {cart['z']}")
        file.write('\n')

