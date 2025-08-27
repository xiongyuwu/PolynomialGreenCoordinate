# PolynomialGreenCoordinate
Generates mesh deformation using cage-based or variational deformation methods based on Polynomial Green Coordinates.
The repository contains the source code for the research paper: http://staff.ustc.edu.cn/~fuxm/projects/CurvedGreen/Poly3DGreen.pdf

## External Libraries

* [Eigen](http://eigen.tuxfamily.org) (Recommended version: 3.3.8)
* [OpenMesh](https://www.openmesh.org) (Recommended version: 10.0 static)
* [Qt](https://www.qt.io) (Recommended version: 5.14.2)
* [Boost](https://www.boost.org) (Recommended version: 1.86.0)

## Compilation

```
git clone -b master git@github.com:xiongyuwu\PolynomialGreenCoordinate.git
cd PolynomialGreenCoordinate
```

Edit lines 14-16 of CmakeLists.txt to set the values of **EIGEN_PATH**,**OPENMESH_PATH** and **OPENMESH_LIB_PATH**.

Edit line 23 of CmakeList.txt to set the value of **CMAKE_PREFIX_PATH** to help find Qt.

Edit line 34 of CmakeList.txt to set the value of **BOOST_ROOT** to help find Boost.

To avoid modifying CMakeLists.txt, you can install the required libraries in the following default locations:

Eigen 3.3.8: D:/ProgramFiles/eigen-3.3.8

OpenMesh 10.0 (Static): D:/ProgramFiles/OpenMesh 10.0 static

Qt 5.14.2: D:/ProgramFiles/Qt/5.14.2

Boost 1.86.0: D:/ProgramFiles/boost_1_86_0

```
mkdir build && cd build
cmake -A x64 ..
```
To view more information and outputs, enable Debug mode:
```
mkdir build && cd build
cmake -A x64 -DDEBUG_MODE=ON ..
```
Open **PolynomialGreenCoordinate.sln**, select **PolynomialGreenCoordinate** as the startup project, and run.

## Usage
Store the cage mesh as "cage.obj" in the \PolynomialGreenCoordinate\build directory.

The default deformed cage uses 3rd-order Bézier facets. To change the facet order, modify the variable "control_point_per_edge" in \PolynomialGreenCoordinate\src\MeshViewer\PGC.h.

For variational deformation, place all constraint files ("fixed_vertices.txt", "position_constraints.txt", "polyline.off", "orientation_constraints.txt") in the \PolynomialGreenCoordinate\build directory.

To generate cage-based deformation using precomputed control points, store the control points data in "deformed_control_points.obj" and place it in the \PolynomialGreenCoordinate\build directory.

- Select File → Open to load the controlled mesh (.obj/.off/.ply/.stl)
- Click **"Change Cage Mode Tri/Quad" button** if the cage is a triangular mesh
- Click **"Read Cage From File" button** to load the cage
- Click **"Draw Control Point" button** to visualize control points
- Click **"Draw Cage" button** to visualize the cage
- Click **"Draw Move Vertices Action" button** to visualize position constraints
- Click **"Generate Deformed Mesh"** for cage-based deformation using current control points
- Click **"Record Origin Mesh" button** to compute and save integral parameters for cage-based deformation
- Click **"Compute Para For Control Points" button** to compute integral parameters for variational deformation
- Click **"Read Integral Para From Files" button** to load precomputed parameters
- Click **"Generate Deformation By Control Points" button** to generate cage-based deformation from "deformed_control_points.obj"
- Click **"Optimize Control Points" button** to generate variational deformation from constraints
- Click **"Generate Variational Deformation" button** to generate deformation with SVD decomposition
- Click **"Generate Para Complete Flow" button** for complete variational deformation workflow
- Click **"Move Control Point" button** for real-time cage deformation by moving control points
- Click **"Move Position Constraints" button** for real-time variational deformation
- Click **"Select Move Control Points" button** to select multiple control points, then use **"Move Control Points" button** to deform
- Click **"Clear Selected Points" button** to clear selections

### Example Flows
**Real-time cage-based deformation (triangular cage):**
1. Place "cage.obj" in \PolynomialGreenCoordinate\build
2. Load controlled mesh
3. Click **"Change Cage Mode Tri/Quad"** 
4. Click **"Read Cage From File"**
5. Click **"Record Origin Mesh"**
6. Click **"Read Integral Para From Files"**
7. Click **"Draw Control Point"** and **"Draw Cage"**
8. Click **"Move Control Point"** and drag control points

**Real-time variational deformation (triangular cage):**
1. Place "cage.obj" and constraint files in \PolynomialGreenCoordinate\build
2. Load controlled mesh
3. Click **"Change Cage Mode Tri/Quad"** 
4. Click **"Read Cage From File"**
5. Click **"Generate Para Complete Flow"**
6. Click **"Draw Move Vertices Action"**
7. Click **"Move Position Constraints"** and drag constraints

## Data Formats
**fixed_vertices.txt**
```
#One line of integers, strore mesh vertex indices that need to fixed
0 1 2 3 4 5 6
```
**position_constraints.txt**

```
#Lines of integers and floats. For each line, the first number is an integer that represents a mesh vertex index, and the following three floats represents the corresponding target position.
0 0.0 0.0 0.0
1 1.0 2.0 3.0
10 3.0 2.0 1.0
```

**polyline.off**

```
#An off file that contains points that represents the medial axis of mesh. This will be used as rigidity constraints in variational deformation.
OFF
0 0 0
0 0 0.5
0 0 1
0 0 1.5
0 0 2
```

**orientation_constraints.txt**

```
#Lines of integers and floats. For each line, the first number is an integer that represents a mesh vertex index, and the following nine floats represents the corresponding target orientaion matrix.
0 1.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 3.0
10 2.0 0.0 0.0 0.0 0.5 -0.866 0 0.866 0.5
```

## Data Sets
Sample datasets are available in the model directory. Each dataset includes a mesh and its control cage. Some also contain deformed control points (for cage-based deformation) or constraints (for variational deformation). These resources enable reproduction of our paper's results.

## Demo
The demo directory contains giraffe.png (reference figure from our paper) and a script to generate its corresponding mesh. Run demo.exe to verify results, with outputs saved in deformed_mesh.obj.

 