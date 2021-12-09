# Reliable Feature-Line Driven Quad-Remeshing

[Nico Pietroni](www.nicopietroni.com), Stefano Nuvoli, 
[Thomas Alderighi](http://vcg.isti.cnr.it/~alderighi/), [Paolo Cignoni](http://vcg.isti.cnr.it/~cignoni/), [Marco Tarini](https://tarini.di.unimi.it/)<br/>
*SIGGRAPH 2021*<br/>

![alt text](teaser.jpg)

## Abstract
We present a new algorithm for the semi-regular quadrangulation of an input surface, driven by its line features, such as sharp creases. We define a perfectly feature-aligned cross-field and a coarse layout of polygonal-shaped patches where we strictly ensure that all the feature-lines are represented as patch boundaries. To be able to consistently do so, we allow non-quadrilateral patches and T-junctions in the layout; the key is the ability to constrain the layout so that it still admits a globally consistent, T-junction-free, and pure-quad internal tessellation of its patches. This requires the insertion of additional irregular-vertices inside patches, but the regularity of the final- mesh is safeguarded by optimizing for both their number and for their reciprocal alignment. In total, our method guarantees the reproduction of feature-lines by construction, while still producing good quality, isometric, pure-quad, conforming meshes, making it an ideal candidate for CAD models. Moreover, the method is fully automatic, requiring no user intervention, and remarkably reliable, requiring little assumptions on the input mesh, as we demonstrate by batch processing the entire Thingi10K repository, with less than 0.5% of the attempted cases failing to produce a usable mesh.

Website: [https://www.quadmesh.cloud/](https://www.quadmesh.cloud/)<br />
DOI: [10.1145/3450626.3459941](https://doi.org/10.1145/3450626.3459941) ACM Transactions on Graphics (SIGGRAPH), 2021

**BibTex**
```
@article{10.1145/3450626.3459941,
author = {Pietroni, Nico and Nuvoli, Stefano and Alderighi, Thomas and Cignoni, Paolo and Tarini, Marco},
title = {Reliable Feature-Line Driven Quad-Remeshing},
year = {2021},
issue_date = {August 2021},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {40},
number = {4},
issn = {0730-0301},
url = {https://doi.org/10.1145/3450626.3459941},
doi = {10.1145/3450626.3459941},
abstract = {We present a new algorithm for the semi-regular quadrangulation of an input surface, driven by its line features, such as sharp creases. We define a perfectly feature-aligned cross-field and a coarse layout of polygonal-shaped patches where we strictly ensure that all the feature-lines are represented as patch boundaries. To be able to consistently do so, we allow non-quadrilateral patches and T-junctions in the layout; the key is the ability to constrain the layout so that it still admits a globally consistent, T-junction-free, and pure-quad internal tessellation of its patches. This requires the insertion of additional irregular-vertices inside patches, but the regularity of the final-mesh is safeguarded by optimizing for both their number and for their reciprocal alignment. In total, our method guarantees the reproduction of feature-lines by construction, while still producing good quality, isometric, pure-quad, conforming meshes, making it an ideal candidate for CAD models. Moreover, the method is fully automatic, requiring no user intervention, and remarkably reliable, requiring little assumptions on the input mesh, as we demonstrate by batch processing the entire Thingi10K repository, with less than 0.5% of the attempted cases failing to produce a usable mesh.},
journal = {ACM Trans. Graph.},
month = {jul},
articleno = {155},
numpages = {17},
keywords = {geometry processing, quad-meshing, modelling}
}
```

### Download
```bash
git clone --recursive https://github.com/nicopietroni/quadwild
```

### Build
Install the libraries boost and gurobi. 
In Ubuntu you can install boost easily with the following terminal commands:
```
apt-get install libboost-dev
```
Open the file libs/libs.pri and set the paths of the requested libraries and the gurobi parameters:
```
#External libraries
BOOST_PATH          = /usr/include/boost/
GUROBI_PATH         = /opt/gurobi950/linux64/
GUROBI_COMPILER     = gurobi_g++5.2
GUROBI_LIB          = gurobi95
```

<br/>

If you do not need CoMISo, you can simply remove the define `COMISO_FIELD` in the file libs.pri:
```
#DEFINES += COMISO_FIELD
```
However, for organic meshes, we suggest to abilitate CoMISo. You need to compile it along with its dependencies (BLAS):
```
apt install libblas-dev
cd quadwild/libs/CoMISo
mkdir build
cd build
cmake ..
make
```

<br/>

You can now compile the project quadwild/quadwild.pro with qmake or QtCreator.

In case you have technical issues or building problems, please write to [stefano.nuvoli@gmail.com](mailto:stefano.nuvoli@gmail.com) or [nico.pietroni@uts.edu.au](mailto:nico.pietroni@uts.edu.au).

### Run
The package is composed of the main command-line quad-remesher (quadwild) and three different components (field_computation, field_tracing, quad_from_patches) that perform different steps of the pipeline.

---

#### quadwild
This project has no visual interface and can be used via command-line. This can be helpful to batch run entire datasets of models. To run the project, once builded, execute the following terminal command:
```
./quadwild <mesh> [.txt setup file] [.rosy file] [.sharp file]
```
The command takes as input a mesh and three optional configuration files:

- **`<mesh>`**: filename of the input triangle mesh. **The mesh can be either an obj or a ply.**
   
- **`.txt setup file` (optional):** The txt setup file contains the parameters in the pipeline. By default, the executable loads the file basic_setup.txt and two other examples are included: basic_setup_mechanical.txt and basic_setup_organic.txt. Any setup parameter can be specified to control the output result. The setup file has the following fields:
```
do_remesh 1 		  //remesh (1) or not (0) the input mesh
sharp_feature_thr 35      //the dihedral angle of sharp features (-1 no features)
alpha 0.02                //regularity vs isometry of the final tessellation. Close to zero -> more regular, Close to 1 -> more singularity are inserted
scaleFact 1               //the scale of the final quadrangulation (the bigger the bigger the quads)
```

- **`.rosy file` (optional)**: This optional file contains parameters for the field computation of the field.
```
fn              //number of faces of the mesh
4               //directions of the field (always 4 for a cross-field)
x0 y0 z0        //XYZ directions of one vector of the cross-field of the first face
...
xn yn zn        //XYZ directions of one vector of the cross-field of the n-th face
```

- **`.sharp file` (optional)**: This optional file contains the informations of the sharp features. Note that, in the pipeline, border edges are considered sharp features by default.
```
sn                //number of sharp features
t0 f0 e0          //for each sharp edge: the first integer is 0 if the edge is concave 1 if convex, then the face and the index of the sharp edge
...
tn fn en          //nth sharp edge
```

The output of quadwild consists of several files:
- **The output smooth quadrangulation (suffix quadrangulation_smooth.obj).**
- The output quadrangulation before being smoothed (suffix quadrangulation.obj).
- Other files:
  - The re-meshed triangulated mesh (suffix rem.obj), the relative field and the sharp features automatically computed (.rosy and .sharp files as above).
  - The mesh decomposed after the tracing (suffix rem_p0.obj).
  - The patch decomposition (.patch file) contains the patch index for each triangle of the rem_p0 mesh.
  - The files .corners, .c_feature, .feature files that contain per patch information (respectively corners of each patch, corners to be fixed and feature lines on the patches).

---

#### field_computation. 
The program can be used either with a GUI or by command line (useful to batch run entire datasets of models).
```
./field_computation <mesh> [.txt setup file] [.rosy file][.sharp file] [batch]
```
The "batch" option makes the program run in the shell without the GUI. The setup file includes additional parameters. By default, the executable loads the file basic_setup.txt.

---

#### field_tracing
This program is used to trace fields and split the mesh into patches.
```
./field_tracing <mesh> [.txt setup file] [batch]
```
It requires having a .rosy and a .sharp file (with the same name of the mesh file). The "batch" option makes the program run in the shell without the GUI. The setup file includes additional parameters. By default, the executable loads the file basic_setup.txt.

---

#### quad_from_patches
This program is used to obtain a quadrangulation from a patch decomposition.
```
./quad_from_patches <mesh> [.txt setup file]
```
It requires to have in the same folder a .corners, .c_feature, .feature files (with the same name of the mesh file). The setup file includes additional parameters. By default, the executable loads the file basic_setup.txt.

## Note
The code has slightly changed and the results could be different from the ones showed in the paper.

## License
[GPL3](LICENSE) licensed
([FAQ](https://www.gnu.org/licenses/gpl-faq.html))



