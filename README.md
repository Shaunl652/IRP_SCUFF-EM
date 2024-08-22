# README

The code in this repository is designed to numerically calculate the information scatter from a dielectric particle trapped in an optical trap. This code is used in the paper [Insert citation for our paper here] and can be used to reproduce the graphs there. We use the scattering solver SCUFF-EM: http://homerreid.github.io/scuff-em-documentation/ to find the scattered file
By providing a new <code>.stl</code> file, the user is also able to calculate the Information Radiation Patterns (IRPs) for any arbitrarily shaped particle. A <code>Makefile</code> has been included to make running the examples easier for the user.

## Files Present

Here we briefly introduce the files in the repository so that the user can make use of them.

1) Various <code>.stl</code> files: These contain the surface meshes for various shaped particles. However, the scattering solver, SCUFF-EM, cannot make use of these. We use a combination of gmsh and mmgs to convert them into <code>.o.msh</code> files that SCUFF-EM can work with. See the mesh commands in the <code>Makefile</code> to see this in action.
2) Various <code>.scuffgeo</code> files: These contain information about the particles in a form that SCUFF-EM can read.
3) <code>Lebedev/</code> folder: This contains various text files with weightings for the Lebedev quadrature rule (See https://en.wikipedia.org/wiki/Lebedev_quadrature) for spherical integration. These values are taken from: https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html and used under the GNU LGPL license.
4) <code>ASM.py</code>: This calculates the Gaussian beam using the Angular Spectrum Method (See Novotny and Hetch PNO Sections 3.4 to 3.6) and makes use of the Lebedev quadrature rule to calculate the spherical integrals.
5) <code>sources_ASM.list</code>: Output from <code>ASM.py</code>, this file contains a list of plane waves that approximate the Gaussian.
6) <code>wavelengths.txt</code>: This contains a list of wavelengths (in microns) at which we want to calculate the IRPs. Currently this only contains the 1.55 microns (1550nm).
7) <code>Sphere.trans</code> and <code>Hex.trans</code>: These files contain the transformations that will be made to the particle in the SCUFF-EM calculations. <code>Sphere.trans</code> only contains translations, while <code>Hex.trans</code> contains the translations and rotations.
8) <code>Eval_Write.py</code>: Writes the file <code>EvalPoints.dat</code> containing the Cartesian coordinates where we measure the scattered field.
9) <code>Fisher_Info.py</code>: Reads in the scattered fields and calculates the IRP based on the Fisher information approach.
10) <code>IRP_laser.py</code>: Reads in the scattered fields and calculates the IRPs using a single incident beam as a reference field.
11) <code>scale.py</code>: Re-scales a <code>.stl</code> file by the scale factor given when calling the code. In the <code>Makefile</code> it is used to re-scale <code>UnitSphere.stl</code> to the radius the user inputs. The returned <code>Sphere.stl</code> is the resized sphere that is then converted to a <code>.o.msh</code> file that SCUFF-EM can then make use of.
12) <code>Hex_Gmsh.py</code>: This can be used to make the <code>.stl</code> files for the hexagonal plates. Simply set the vertex-to-vertex diameter and thickness by changing the variables <code>D</code> and <code>T</code> respectively.

## Example

Here we give a brief example for finding the IRP of a sphere.

1) In your command line run the command <code>Eval_Points</code>. This will make the <code>EvalPoints.dat</code> containing the locations where the field should be calculated.
2) Next run the command <code>make mesh_Sphere</code>.
3) You will now be asked to input the radius of the sphere in microns. To achieve the IRP at $R/\lambda = 0.01$ the radius you input would be <code>0.0155</code>.
4) You will now be asked for the Hausdorff distance. This is the maximum distance from the perfect sphere to the re-meshed version (see https://en.wikipedia.org/wiki/Hausdorff_distance). I.e. a smaller value will lead to a more dense mesh. When the re-meshing is compleate, Gmsh will show the final mesh and you can go back to step 1) to make adjustments as necessary.
5) When you are happy with the mesh, run the command <code>Simulate_Sphere</code>. This will run SCUFF-EM with the relevant args file. This will produce a <code>.scattered</code> and a <code>.total</code> file which contains the scattered and total fields at each point in <code>EvalPoints.dat</code>.
6) Once you have the <code>.scattered</code> file, you can plot the IRPs. To do this make sure the <code>File_Name</code> variable in <code>Fisher_Info.py</code> is correct. This should be the first part of the <code>.scattered</code> file. In this case the file will be <code>Sphere.EvalPoints.scattered</code> so in <code>Fisher_Info.py</code> we will set <code> File_Name = 'Sphere'</code>. Also ensure the variable <code>Rotation</code> is correct. If True, the code will attempt to find the IRPs for the rotational degrees of freedom. For this case an error will be returned, but that translational DoFs will have been calculated already, so this is not a problem.
7) You can now find your IRPs in the <code>Plots/</code> folder.
