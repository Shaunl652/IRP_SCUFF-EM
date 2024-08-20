import trimesh
from sys import argv
f1 = argv[1]
s = float(argv[2])
f2 = argv[3]

m = trimesh.load(f1)
m.vertices *= s
m.export(f2)
