# Scotty3D
Scotty3D for mesh edit

Mesh edit
KEY COMMANDS:
Key					Action
space				switch to navigate mode
tab					show/hide info panel
hold control		temporarily switch to translation in edit mode
hold shift			temporarily switch to scaling in edit mode
hold alt/option		temporarily switch to rotation in edit mode
e					cycle through (e)dit modes
b					toggle bevel mode
n					select next halfedge
t					select twin halfedge
h					select halfedge of current element
T					Triangulate mesh
s					subdivide Catmull-Clark
S					Subdivide linear
backspace/delete	erase selected edge
f					flip selected edge
c					collapse selected edge
p					split selected edge triangle meshes only!
u					upsample triangle meshes only!
i					isotropic remesh triangle meshes only!
d					downsample triangle meshes only!
w 					then 0--9	write scene to numbered buffer
l 					then 0--9	load scene from numbered buffer**



PathTracer
You can run pathtracer with command:  ./scotty3d -s 16 -m 4 -t 8 ../../dae/sky/CBspheres.dae on Windows.
KEY COMMANDS:
Command										Key
Return to mesh edit mode					M
Show BVH visualizer mode					V
Show ray traced output						R
Decrease area light samples (RT mode)		-
Increase area light samples (RT mode)		+
Decrease samples (camera rays) per pixel	[
Increase samples (camera rays) per pixel	]
Descend to left child (BVH viz mode)		<
Descend to right child (BVH viz mode)		>
Move to parent node (BVH viz mode)			?
Reset camera to default position			SPACE
Edit a vertex position						(left-click and drag on vertex)
Rotate camera								(left-click and drag on background)
Zoom camera									(mouse wheel)
Dolly (translate) camera					(right-click and drag on background)


