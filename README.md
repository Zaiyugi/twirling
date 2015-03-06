Horn Modeling
==============
Created: 2014-04-30
Updated: 2015-02-11

Project Description
====================
3D Modeling tool implementing and exploring techniques put forth by the paper Twirling Sculptures by Ergun Akleman. Built using OpenMesh and OpenGL

Building
=========

| Command         | Use |
| --------------- | --- |
| **make all**    | Build and create the lib files |
| **make launch** | Build the main code |
| **./launch**    | Run |

Note on building:
OpenMesh is used for writing and storing meshes within the program. Depending on how it is installed, you might need to modify the OPENMESH variables in the Makefile to have the right paths. Everything else is sourced locally.

During operation
=================

* Subdivisions control surface subdivisions; Subdivided Bezier not included. Doo-Sabin subdivision is used and should work for all cases

* When Importing/Exporting, the Surface Name field is the path/name of the file, minus the extension. So "models/icosa.obj" becomes "models/icosa". 

* Extrusion is controlled by selection and the type. Each type will consider only the faces selected, so to extrude all faces, select all faces. Otherwise, switch to 'Select' mode and select the desired region. Extrusion runs per face, not per region; extruding a region of faces is no different than extruding each face individually. 

* When extruding along a path, then step size and step count parameters are ignored; the number of steps and size of each is determined by the number of points and distance between them on the path.

* I added several editing commands; they're grouped under main controls. They should make manipulating surfaces easier. A history is also kept, so you can undo some commands. 

* Both extrusion and subdivision rely on application. To subdivide, set the number of subdivisions desired and press 'Apply Subdivisions'. A new mesh, reflecting that subdivision, will be pushed onto the top of the history. The same applies to extrusion.
