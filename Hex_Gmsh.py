#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
#
#  Building a hexagonal disk from user entered values
#
#  Based on tutorial 1 on the Gmsh website: https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_11_1/tutorials/python/t1.py
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys
from numpy import sqrt

D = 5      # Winstone et al. gives plate diameter in microns, use this for ease
R = D/2    # Radius in microns (distance from center to vertex)
T = 0.2/2  # disk thickness in microns (distance from centre to each face)



# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# Next we add a new model named "Hex" 
gmsh.model.add("Hex")


# We first create the points of the verticies
# - the first 3 arguments are the point coordinates (x, y, z)
# - the next (optional) argument is the target mesh size close to the point
# - the last (optional) argument is the point tag (a stricly positive integer
#   that uniquely identifies the point)
lc = 1e-2
# Now we just want to define all the points on the front face
B = gmsh.model.geo.addPoint( R      , 0          , T, lc)
C = gmsh.model.geo.addPoint( R-(R/2), R*sqrt(3)/2, T, lc)
D = gmsh.model.geo.addPoint(-R+(R/2), R*sqrt(3)/2, T, lc)
E = gmsh.model.geo.addPoint(-R      , 0          , T, lc)
F = gmsh.model.geo.addPoint(-R+(R/2),-R*sqrt(3)/2, T, lc)
G = gmsh.model.geo.addPoint( R-(R/2),-R*sqrt(3)/2, T, lc)


# Next we do the same thing for the back face
H = gmsh.model.geo.addPoint( R      , 0          ,-T, lc)
I = gmsh.model.geo.addPoint( R-(R/2), R*sqrt(3)/2,-T, lc)
J = gmsh.model.geo.addPoint(-R+(R/2), R*sqrt(3)/2,-T, lc)
K = gmsh.model.geo.addPoint(-R      , 0          ,-T, lc)
L = gmsh.model.geo.addPoint(-R+(R/2),-R*sqrt(3)/2,-T, lc)
M = gmsh.model.geo.addPoint( R-(R/2),-R*sqrt(3)/2,-T, lc)

# Now we have to define our edges by building lines bwtween the appropriate vertices
# Start with the front face
BC = gmsh.model.geo.addLine(B,C)
CD = gmsh.model.geo.addLine(C,D)
DE = gmsh.model.geo.addLine(D,E)
EF = gmsh.model.geo.addLine(E,F)
FG = gmsh.model.geo.addLine(F,G)
GB = gmsh.model.geo.addLine(G,B)


# And again for the back face
HI = gmsh.model.geo.addLine(H,I)
IJ = gmsh.model.geo.addLine(I,J)
JK = gmsh.model.geo.addLine(J,K)
KL = gmsh.model.geo.addLine(K,L)
LM = gmsh.model.geo.addLine(L,M)
MH = gmsh.model.geo.addLine(M,H)

# Now we join the faces
BH = gmsh.model.geo.addLine(B,H)
CI = gmsh.model.geo.addLine(C,I)
DJ = gmsh.model.geo.addLine(D,J)
EK = gmsh.model.geo.addLine(E,K)
FL = gmsh.model.geo.addLine(F,L)
GM = gmsh.model.geo.addLine(G,M)

# Before we can define the surfaces, we must define a set of loops to contain the surface
loop_front = gmsh.model.geo.addCurveLoop([BC,CD,DE,EF,FG,GB])
loop_back  = gmsh.model.geo.addCurveLoop([HI,IJ,JK,KL,LM,MH])

# Can't forget about the surfaces on the sides
loop_side1 = gmsh.model.geo.addCurveLoop([ BC, CI,-HI,-BH])
loop_side2 = gmsh.model.geo.addCurveLoop([ CD, DJ,-IJ,-CI])
loop_side3 = gmsh.model.geo.addCurveLoop([ DE, EK,-JK,-DJ])
loop_side4 = gmsh.model.geo.addCurveLoop([ EF, FL,-KL,-EK])
loop_side5 = gmsh.model.geo.addCurveLoop([ FG, GM,-LM,-FL])
loop_side6 = gmsh.model.geo.addCurveLoop([ GB, BH,-MH,-GM])


# We can then define the surface as a list of curve loops
surf_front = gmsh.model.geo.addPlaneSurface([loop_front])
surf_back  = gmsh.model.geo.addPlaneSurface([loop_back])

# And the sides
surf_side1 = gmsh.model.geo.addPlaneSurface([loop_side1])
surf_side2 = gmsh.model.geo.addPlaneSurface([loop_side2])
surf_side3 = gmsh.model.geo.addPlaneSurface([loop_side3])
surf_side4 = gmsh.model.geo.addPlaneSurface([loop_side4])
surf_side5 = gmsh.model.geo.addPlaneSurface([loop_side5])
surf_side6 = gmsh.model.geo.addPlaneSurface([loop_side6])

# Before they can be meshed the surfaces must be syncrinised
gmsh.model.geo.synchronize()


# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# ... and save it to disk
gmsh.write("Hex.stl")



# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
# command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
