# Scotty3D

## Scotty3D for Mesh Edit

### Mesh Edit

#### Key Commands:

| Key              | Action                                         |
| ---------------- | ---------------------------------------------- |
| space            | Switch to navigate mode                        |
| tab              | Show/hide info panel                           |
| hold control     | Temporarily switch to translation in edit mode |
| hold shift       | Temporarily switch to scaling in edit mode     |
| hold alt/option  | Temporarily switch to rotation in edit mode    |
| e                | Cycle through (e)dit modes                     |
| b                | Toggle bevel mode                              |
| n                | Select next halfedge                           |
| t                | Select twin halfedge                           |
| h                | Select halfedge of current element             |
| T                | Triangulate mesh                               |
| s                | Subdivide Catmull-Clark                        |
| S                | Subdivide linear                               |
| backspace/delete | Erase selected edge                            |
| f                | Flip selected edge                             |
| c                | Collapse selected edge                         |
| p                | Split selected edge (triangle meshes only!)    |
| u                | Upsample (triangle meshes only!)               |
| i                | Isotropic remesh (triangle meshes only!)       |
| d                | Downsample (triangle meshes only!)             |
| w then 0--9      | Write scene to numbered buffer                 |
| l then 0--9      | Load scene from numbered buffer                |

---

## PathTracer

You can run the path tracer with the command: ./scotty3d -s 16 -m 4 -t 8 ../../dae/sky/CBspheres.dae

This works on Windows.

#### Key Commands:

| Command                                  | Key                                |
| ---------------------------------------- | ---------------------------------- |
| Return to mesh edit mode                 | M                                  |
| Show BVH visualizer mode                 | V                                  |
| Show ray traced output                   | R                                  |
| Decrease area light samples (RT mode)    | -                                  |
| Increase area light samples (RT mode)    | +                                  |
| Decrease samples (camera rays) per pixel | [                                  |
| Increase samples (camera rays) per pixel | ]                                  |
| Descend to left child (BVH viz mode)     | <                                  |
| Descend to right child (BVH viz mode)    | >                                  |
| Move to parent node (BVH viz mode)       | ?                                  |
| Reset camera to default position         | SPACE                              |
| Edit a vertex position                   | Left-click and drag on vertex      |
| Rotate camera                            | Left-click and drag on background  |
| Zoom camera                              | Mouse wheel                        |
| Dolly (translate) camera                 | Right-click and drag on background |

---
