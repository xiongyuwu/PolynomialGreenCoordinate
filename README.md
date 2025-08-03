#PolynomialGreenCoordinate
Generate mesh deformation with cage-based deformation or variational deformation based on polynomial green coordinates. 


The program is used to load meshes and render using OpenGL.

## External Libraries

* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/), Recommended version: the latest 8.1(at Oct. 2020)
* [Qt](https://www.qt.io/), Recommended version: 5.13.0

## Compilation

```
git clone -b master git@github.com:xiongyuwu/PolynomialGreenCoordinate.git
cd PolynomialGreenCoordinate
```

Edit lines 20-23 of CmakeLists.txt to set the values of **Qt5Core_DIR**,**Qt5Gui_DIR**,**Qt5OpenGL_DIR** and **Qt5Widgets_DIR**
Edit lines 31-33 of CmakeLists.txt to set the values of **EIGEN_PATH**,**OPENMESH_PATH** and **OPENMESH_LIB_PATH**

```
mkdir build && cd build
cmake -A x64 ..
```
use below code to show more information
```
mkdir build && cd build
cmake -A x64 -DDebugMode=ON ..
```

Open **PolynomialGreenCoordinate.sln**, select **PolynomialGreenCoordinate** as launch project, and run.

##Usage
Store cage mesh as "cage.obj", then put it in /PolynomialGreenCoordinate/build directory.
If you want to run variational deformation, also put all constraints("fixed_vertices.txt", "position_constraints.txt", "polyline.off", "orientation_constraints.txt") in /PolynomialGreenCoordinate/build directory.
If you want to generate cage-based deformation by reading deformed control points, store the control points data in "deformed_control_points.obj", then put it in /PolynomialGreenCoordinate/build directory.

Select File->Open to and load controlled mesh (.obj .off .ply .stl).
Click "Change Cage Mode Tri/Quad" buttom if the cage mesh is a triangular mesh.
Click "Read Cage From File" buttom to read cage.
Click "Draw Control Point" buttom to visualize control points.
Click "Draw Cage" buttom to visualize cage.
Click "Draw Move Vertices Action" buttom to visualize position constraints.
Click "Generate Deformed Mesh" to generate cage-based deformation based on current control points.
Click "Record Origin Mesh" buttom to compute and write integral paras for cage-based deformation.
Click "Compute Para For Control Points" buttom to compute integral paras of mesh vertices for variational deformation.
Click "Read Integral Para From Files" buttom to read precompted integral paras in /PolynomialGreenCoordinate/build directory.
Click "Generate Deformation By Control Points" buttom to read "deformed_control_points.obj" and generate corresponding cage-based deformation.
Click "Optimize Control Points" buttom to read constraints and precompute integrals and generate variational deformation.
Click "Generate Variational Deformation" buttom to read constraints, precompute integrals and precompute SVD decomposition and generate variational deformation.
Click "Generate Para Complete Flow" buttom to generate complete flow of variational deformation, including reading constarints, computing integral paras and generate deformed mesh.
Click "Move Control Point" buttom to generate real-time cage-based deformation by moving control points. 
Click "Move Position Constraints" buttom to generate real-time variational deformation by moving position constraints. 
Click "Select Move Control Points" buttom to select multiple control points, then click "Move Control Points" buttom to generate cage-based deformation while moving these control points together.
Click "Clear Selected Points" buttom to clear selected control points.

#Example flow
real-time cage-based deformation with triangular cage:
Put "cage.obj" in /PolynomialGreenCoordinate/build directory->load controlled mesh->click "Change Cage Mode Tri/Quad" buttom to use triangular cage->click "Record Origin Mesh" buttom to compute and write integral paras->click "Read Integral Para From Files" buttom to load precomputed integral datas->Click "Draw Control Point" and "Draw Cage" buttom to visualize cage and control points->Click "Move Control Point" buttom and drag control points to generate real-time cage-based deformation.

real-time variational deformation with triangular cage:
Put "cage.obj", "fixed_vertices.txt", "position_constraints.txt", "polyline.off", "orientation_constraints.txt" in /PolynomialGreenCoordinate/build directory->load controlled mesh->click "Change Cage Mode Tri/Quad" buttom to use triangular cage->click "Generate Para Complete Flow" buttom to generate variational deformation based on input constarints->click "Draw Move Vertices Action" buttom to visualize position constraints->click "Move Position Constraints" buttom and drag position constraints to generate real-time variational deformation.

##Data Format
Here are required format of input constraints.
```fixed_vertices.txt
#One line of integers, strore mesh vertex indices that need to fixed
0 1 2 3 4 5 6
```
```position_constraints.txt
#Lines of integers and floats. For each line, the first number is an integer that represents a mesh vertex index, and the following three floats represents the corresponding target position.
0 0.0 0.0 0.0
1 1.0 2.0 3.0
10 3.0 2.0 1.0
```
```polyline.off
#An off file that contains points that represents the medial axis of mesh. This will be used as rigidity constraints in variational deformation.
OFF
0 0 0
0 0 0.5
0 0 1
0 0 1.5
0 0 2
```
```orientation_constraints.txt
#Lines of integers and floats. For each line, the first number is an integer that represents a mesh vertex index, and the following nine floats represents the corresponding target orientaion matrix.
0 1.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 3.0
10 2.0 0.0 0.0 0.0 0.5 -0.866 0 0.866 0.5
```

