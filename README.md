# README

The code in this repository is designed to numerically calculate the information scatter from a dielectric particle trapped in an optical trap. This code is used in the paper [Insert citation for our paper here] and can be used to reproduce the graphs there. We use the scattering solver SCUFF-EM: http://homerreid.github.io/scuff-em-documentation/ to find the scattered file
By providing a new <code>.stl</code> file, the user is also able to calculate the Information Radiation Patterns (IRPs) for any arbitrarily shaped particle. A <code>Makefile</code> has been included to make running the examples easier for the user.

## Files Present

Here we briefly introduce the files in the repository so that the user can make use of them.

1) Various <code>.stl</code> files: These contain the surface meshes for various shaped particles. However, the scattering solver, SCUFF-EM, cannot make use of these. We use a combination of gmsh and mmgs to convert them into <code>.o.msh</code> files that SCUFF-EM can work with. See the mesh commands in the <code>Makefile</code> to see this in action.
2) Various <code>.scuffgeo</code> files: These contain information about the particles in a form that SCUFF-EM can read.
