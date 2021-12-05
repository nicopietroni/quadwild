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
git clone --recursive https://github.com/stefanonuvoli/quadwild
```

### Build
Install the libraries boost and gurobi. 
In Ubuntu you can install boost easily with the following terminal commands:
```
sudo apt-get install libboost-dev
```
Open the file libs/libs.pri and set the paths of the requested libraries and the gurobi parameters:
```
#External libraries
BOOST_PATH          = /usr/include/boost/
GUROBI_PATH         = /opt/gurobi903/linux64/
GUROBI_COMPILER     = gurobi_g++5.2
GUROBI_LIB          = gurobi90
```
You can now compile the project quadwild/quadwild.pro with qmake or QtCreator.
<br /><br/>
In case you have technical issues or building problems, please write to [stefano.nuvoli@gmail.com](mailto:stefano.nuvoli@gmail.com) or [nico.pietroni@uts.edu.au](mailto:nico.pietroni@uts.edu.au).

## Run
The package is composed of the main command line quadrangular, which is in the folder quad wild and the three different components of the pipeline.

### QuadWild
The program has no visual interface and can be used via command-line (this is valuable to batch run entire datasets of models).
```
./quadwild <mesh> [.txt setup file] [.rosy file][.sharp file]
```
#### The mesh can be either an obj or a ply.

#### txt setup file (optional)
By default, the executable load the file basic_setup.txt and two other examples are included, named basic_setup_mechanical.txt and basic_setup_organic.txt. Any setup parameter can be specified to control the output. The setup file has the following fields:
```
do_remesh 1/0 		  // remesh or not the input mesh
sharp_feature_thr 35      // the dihedral angle of sharp features (-1 no features)
alpha 0.02                //regularity vs isometry of the final tessellation. Close to zero -> more regular, Close to 1 -> more singularity are inserted
scaleFact 1               //the scale of the final quadrangulation (the bigger the bigger the quads)
```
#### Rosy file format (optional)
```
fn              // number of faces of the mesh
4               // directions of the field (always 4 for a cross-field)
x0 y0 z0        // XYZ directions of one vector of the cross-field of the first face
...
xn yn zn        // XYZ directions of one vector of the cross-field of the n-th face
```

#### Sharp file format (optional)
```
sn                // number of sharp features
t0 f0 e0          // for each sharp edge: the first integer is 0 if the edge is concave 1 if convex, then the face and the index of the sharp edge
...
tn fn en          // nth sharp edge
```
Notice that border edges are considered sharp features by default

#### Output
The program outputs:

- Two quadrangulation files. One that has been smoothed and one not.

Other:

- The re-meshed triangulated mesh (suffix rem.obj), the relative field and the sharp features (.rosy and .sharp files as above)
- The mesh decomposed after the tracing suffix rem_p0.obj)
- The patch decomposition (.patch file) contains the patch index for each triangle of the rem_p0 mesh.
.corners, .c_feature, .feature files that contain per patch information (corners of each patch, corners to be fixed and feature lines on the patches)

### Components

#### field_computation. The program can be used either with a GUI or by command line (useful to batch run entire datasets of models).

```
./field_computation <mesh> [.txt setup file] [.rosy file][.sharp file] [batch]
```

The "batch" option makes the program run in the shell without the GUI. The setup file includes additional parameters. The one loaded by default is basic_setup.txt.


#### field_tracing
This program is used to trace fields and split the mesh into patches.

```
./field_tracing <mesh> [.txt setup file] [batch]
```
It requires having a .rosy and a .sharp file (with the same name of the mesh file). The "batch" option makes the program run in the shell without the GUI. The setup file includes additional parameters. The one loaded by default is basic_setup.txt.

#### quad_from_patches
This program is used to obtain a quadrangulation from a patch decomposition.

```
./quad_from_patches <mesh> [.txt setup file]
```
It requires to have in the same folder a .corners, .c_feature, .feature files
 (with the same name of the mesh file). The setup file includes additional parameters. The one loaded by default is basic_setup.txt.

The main command
## Note
The code has slightly changed and the results could be different from the ones showed in the paper.

## License
[GPL3](LICENSE) licensed
([FAQ](https://www.gnu.org/licenses/gpl-faq.html))



