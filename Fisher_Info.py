# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:19:17 2023

@author: shaun
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Reads in the scattered field for the various points and the reference field
# Then performs the process as described in the Novotny paper to find the IRP plots
# However, this code does not use a refrence field, instead we follow the method 
# given in https://arxiv.org/pdf/2309.00010.pdf Eq.(4)
# =============================================================================

import numpy as np

from scipy.interpolate import griddata
from scipy.constants import speed_of_light, epsilon_0, hbar, pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from cmath import phase
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path


# These are the varibles to edit

File_Name = 'Sphere' # Name of the data file
Rotation  = True # Set to true if you ant to see the roational info

# Sets up the NAs to calculate the detection efficencies
NA_L = 0.75
NA_R = 0.75
theta_L = np.arcsin(NA_L) # Angle covered by LEFT lens
theta_R = np.arcsin(NA_R) # Angle covered by RIGHT lens

# These first two functions are mostly about reading in the data
def valid(line):
    # Returns False if the first character in the line is '#'
    if line.startswith("#"): return False
    return True

def extract(line):
    # Pulls the values out of the file
    values = line.strip().split()
    return [v for v in values]

def points(r,t,p):
    from numpy import sin, cos
    # Converts the spherical polars to cartesian coords
    x = r*sin(t)*cos(p)
    y = r*sin(t)*sin(p)
    z = r*cos(t)
    
    return x,y,z

def build_E_field(data,tag):
    # Takes in the feild components at each location and converts them to an
    # array of complex number
    output = []
    # Loops through each line in the data set (maybe a quicker way to do this, but I'll find it later (?))
    for line in data:
        # Only does the calculations for the right translation
        if line[4] == tag:
            Ex = float(line[6])  + 1j*float(line[7])
            Ey = float(line[8])  + 1j*float(line[9])
            Ez = float(line[10]) + 1j*float(line[11])
            output.append([Ex,Ey,Ez])
    return np.array(output)#/1e-6

def build_H_field(data,tag):
    # Takes in the feild components at each location and converts them to an
    # array of complex number
    output = []
    # Loops through each line in the data set (maybe a quicker way to do this, but I'll find it later (?))
    for line in data:
        # Only does the calculations for the right translation
        if line[4] == tag:
            Hx = float(line[12]) + 1j*float(line[13])
            Hy = float(line[14]) + 1j*float(line[15])
            Hz = float(line[16]) + 1j*float(line[17])
            output.append([Hx,Hy,Hz])
    return np.array(output)#/1e-6

# Differential approximation
def dfdx(f,x):
    return (f[1]-f[0])/(x[1]-x[0])

# Keep all of these in microns for now
WL = 1550e-9 # Laser wavelength in m
k = 2*pi/WL # Wave number
W0 = 6 # Laser waist in um
E0 = 1 # incident E-field strength 
omega = k*speed_of_light # angluar freeq of laser

# Depends if we want to consider rotational DoFs
if Rotation:
    axis_list = ['x','y','z','thetax','thetay','thetaz']
    
else:
    axis_list = ['x','y','z']

# Finally we read in the data
with open(f'{File_Name}.EvalPoints.scattered','r') as file:
    lines = (line.strip() for line in file if valid(line))
    data_sca = np.array([extract(line) for line in lines])
    
    
# These are the cartesian coordinates of where we measure the fields
coords_str = [vals[:3] for vals in data_sca if vals[4]=='-x']
x = []
y = []
z = []

for point in coords_str:
    x.append(point[0])
    y.append(point[1])
    z.append(point[2])
    
# Converts the measurment locations to floats
xs = [float(xi) for xi in x]
ys = [float(yi) for yi in y]
zs = [float(zi) for zi in z]

# Sets up an array containing the measurment locations
coords = np.array([xs,ys,zs]).T/1e6

# Set up the spherical polars
# These first few lines are just picking equal spacing on a unit sphere (see: https://mathworld.wolfram.com/SpherePointPicking.html)
# We use more points than we measured and hope for interpolation
u = np.linspace(0,1,200)
v = np.linspace(1,0,200)
theta = np.arccos(2*v-1)
phis = 2*np.pi*u


# Start on detection efficeny stuff
# This block of code finds the number of array elements to cover the lens NA

# The difference between each point in the theta space
info_dtheta = np.gradient(theta)

# We only look at theta as we're assuming a circular lens, so we go over theta_NA
# in the theta axis and 2pi in the phi axis

# Finds the number of elements for the integral in the left detection efficency calculation
R_count = 0
angle_R = 0

while angle_R < theta_R:
    angle_R += info_dtheta[R_count]
    R_count+=1
    
# Finds the number of elements for the integral in the right detection efficency calculation
# The LEFT lens is centered on theta = pi which is the end of the array
L_count = len(info_dtheta)-1
angle_L = 0

while angle_L < theta_L:
    angle_L += info_dtheta[L_count]
    L_count-=1
    

# Mesh grid these things to get 2D polts of the angles
THETA, PHI = np.meshgrid(theta, phis)
# Find how x,y, and z vary across the unit sphere
X,Y,Z = points(1,THETA,PHI)

dTHETA = np.gradient(-np.cos(THETA),axis=1)

dPHI = np.gradient(PHI,axis=0)
#Unit_Vec = np.array([np.cos(PHI)*np.sin(THETA),np.sin(PHI)*np.sin(THETA),np.cos(THETA)])
domega = dTHETA*dPHI # Builds the sin(theta part for the integrations over all the spehre)

    
flux_figs = {} # Save the figures so we can look at them in the matplot window
eta_vals = {}  # Save the detection efficencies so we can see them


def Fisher_Info(mu,angle=False):
    # This function finds the Fisher info for the mu DoF at each point
    
    # Get the +ve and -ve translations for each axis
    trans_str = ['-'+mu,'+'+mu]
    
    # If angle is true, then we are considering a rotational DoF
    if angle:
        # Rotate by 1 degree in each direction
        trans_flt = [-1,+1]
    else:
        # Move by lambda/100 in each direction
        trans_flt = [-0.0155e-6,+0.0155e-6]
    
    # Scattered feilds at each detector
    E_Field = {loc: build_E_field(data_sca,loc)  for loc in trans_str}
    H_Field = {loc: build_H_field(data_sca,loc)  for loc in trans_str}
    
    S_FI_temp = [] # stores the FI flux at each detector
    
    for det_pos,x in enumerate(xs):
        print(f'\rDetector {det_pos+1:05}/{len(xs)}',end='')
        Ex_x = []
        Ey_x = []
        Ez_x = []
        Hx_x = []
        Hy_x = []
        Hz_x = []
        # Finds the power at the detector as a function of satterer position
        for trans in trans_str:
            # Finds the total x,y,z components for the field at the dectecor for
            # for the detector at 'det_pos' with the particle at location 'trans'
            Ex,Ey,Ez = np.conj(E_Field[trans][det_pos])
            Ex_x.append(Ex)
            Ey_x.append(Ey)
            Ez_x.append(Ez)
            
            Hx,Hy,Hz = H_Field[trans][det_pos]
            Hx_x.append(Hx)
            Hy_x.append(Hy)
            Hz_x.append(Hz)

        # Usimple derivitive works for linear change in field (small changes in mu)
        dEx = dfdx(Ex_x,trans_flt)
        dEy = dfdx(Ey_x,trans_flt)
        dEz = dfdx(Ez_x,trans_flt)
    
        dHx = dfdx(Hx_x,trans_flt)
        dHy = dfdx(Hy_x,trans_flt)
        dHz = dfdx(Hz_x,trans_flt)
        
        # Puts the components into vectors
        dEdx = np.array([dEx,dEy,dEz])
        dHdx = np.array([dHx,dHy,dHz])

        # makes the fisher info at each detector
        S_FI_vec = 2*np.real(np.cross(dEdx,dHdx))/(hbar*omega)
        Unit_Vec = coords[det_pos]
        S_FI_temp.append(S_FI_vec@Unit_Vec)
        
    return np.array(S_FI_temp)

from scipy import interpolate
u_interp = np.linspace(0,1,500)
v_interp = np.linspace(1,0,500)
theta_interp = np.arccos(2*v-1)
phis_interp = 2*np.pi*u


# Now we work through each DoF
for motionaxis in axis_list:
    print(f'Working on motion in the {motionaxis} axis')
    
    if motionaxis=='thetax' or motionaxis=='thetay' or motionaxis=='thetaz':
        S_FI = Fisher_Info(motionaxis,angle=True)
    else:
        S_FI = Fisher_Info(motionaxis,angle=False)
    # Remeshes the spectral density from a list of locations to a 2D grid at the locations in X,Y,Z
    S_FI_grid = griddata(coords, np.array(S_FI), (X,Y,Z), method='nearest')
    
    f = interpolate.interp2d(theta,phis,S_FI_grid)
    S_FI_interp = f(theta_interp,phis_interp)
    
    # Finds the information radiated to each angle (See the start of Sec II C)
    I = S_FI_interp/np.trapz(np.trapz(S_FI_interp*np.sin(THETA),theta,axis=1),phis) #
    
    print('\nFound Info Plot')
    print('Finding detection efficiencies')
    
    # So, we want to integrate over the angle theta_L
    # int_0^2pi int_0^theta_L I sin(theta) dtheta dphi

    # Detection effiencies for the RIGHT and LEFT lens respectivly
    eta_R = np.trapz(np.trapz((I[:,:R_count+1]*np.sin(THETA[:,:R_count+1])).T,theta[:R_count+1],axis=0),phis)
    
    eta_L = np.trapz(np.trapz((I[:,L_count:]*np.sin(THETA[:,L_count:])).T,theta[L_count:],axis=0),phis)

    # Save the etas
    eta_vals[motionaxis] = {'left': eta_L, 'right': eta_R}
    
    # Put the plots together
    Xi = X*I
    Yi = Y*I
    Zi = Z*I
    d = np.sqrt(Xi**2+Yi**2+Zi**2)
    d_norm = d/d.max()
    
    # Plots the information patteren
    plot_path = f'Plots/{File_Name}'
    Path(plot_path).mkdir(parents=True, exist_ok=True) 
    
    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
    
    ax.plot_surface(Xi,Yi,Zi,facecolors=plt.cm.jet(d_norm),shade=False)
    ax.set(xlabel='x',ylabel='y',zlabel='z')
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(d)
    #fig.colorbar(m,ax=ax)
    
    # Set the view angle
    ax.view_init(elev=20,azim=30,vertical_axis='x')
    
    # Set the axis limmits to be the same to see the right shape
    max_range = np.array([Xi.max()-Xi.min(), Yi.max()-Yi.min(), Zi.max()-Zi.min()]).max() / 2.0
    
    # Finds the mid points
    mid_x = (Xi.max()+Xi.min()) * 0.5
    mid_y = (Yi.max()+Yi.min()) * 0.5
    mid_z = (Zi.max()+Zi.min()) * 0.5
    
    # Set the limits
    Lim_Mod = 0.8
    xmin = (mid_x - max_range)*Lim_Mod
    xmax = (mid_x + max_range)*Lim_Mod
    ymin = (mid_y - max_range)*Lim_Mod
    ymax = (mid_y + max_range)*Lim_Mod
    zmin = (mid_z - max_range)*Lim_Mod
    zmax = (mid_z + max_range)*Lim_Mod
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(f'{plot_path}/{motionaxis}.pdf')
    fig.savefig(f'{plot_path}/{motionaxis}.png')
    flux_figs[motionaxis] = fig
    plt.close()


















