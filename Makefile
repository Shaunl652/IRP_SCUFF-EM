Simulate_Hex: Hex.o.msh
	-rm Hex.EvalPoints.scattered
	-rm Hex.EvalPoints.total
	scuff-scatter < args_Hex

Simulate_Rod: Rod.o.msh
	-rm Rod.EvalPoints.total
	-rm Rod.EvalPoints.scattered
	scuff-scatter < args_Rod

Simulate_Sphere: Sphere.o.msh
	-rm Sphere.EvalPoints.total
	-rm Sphere.EvalPoints.scattered
	scuff-scatter < args_Sphere

mesh_Sphere: UnitSphere.stl
	@read -p "Enter particle radius in microns: " Size; \
	python3 scale.py UnitSphere.stl $$Size Sphere.stl
	# This line sets the correct format
	gmsh Sphere.stl -2 -format msh2
	# This line re-meshes the mesh and controls the mesh density
	@read -p "Enter Hausdorff distance: " Mesh; \
	mmgs -hausd $$Mesh -nr Sphere.msh
	gmsh Sphere.o.msh

mesh_Hex: Hex.stl
	gmsh Hex.stl -2 -format msh2
	@read -p "Enter Hausdorff distance: " Mesh; \
	mmgs -hausd $$Mesh -nr Hex.msh
	gmsh Hex.o.msh

mesh_Rod: Rod.stl
	gmsh Rod.stl -2 -format msh2
	@read -p "Enter Hausdorff distance: " Mesh; \
	mmgs -hausd $$Mesh -nr Rod.msh

Eval_Points:
	# Makes the points at which to measure the EM field
	python3 Eval_Write.py

clean:
	# Cleans out all the old files
	-rm *.msh
	-rm *.log
	-rm *.pft
	-rm *.scattered
	-rm *.total

