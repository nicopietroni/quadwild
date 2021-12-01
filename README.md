# Reliable Feature-Line Driven Quad-Remeshing

[Nico Pietroni](https://profiles.uts.edu.au/Nico.Pietroni), Stefano Nuvoli, 
[Thomas Alderighi](http://vcg.isti.cnr.it/~alderighi/), [Paolo Cignoni](http://vcg.isti.cnr.it/~cignoni/), [Marco Tarini](https://tarini.di.unimi.it/)<br/>
*SIGGRAPH 2021*<br/>

![alt text](teaser.png)

## Abstract
We present a new algorithm for the semi-regular quadrangulation of an input surface, driven by its line features, such as sharp creases. We define a perfectly feature-aligned cross-field and a coarse layout of polygonal-shaped patches where we strictly ensure that all the feature-lines are represented as patch boundaries. To be able to consistently do so, we allow non-quadrilateral patches and T-junctions in the layout; the key is the ability to constrain the layout so that it still admits a globally consistent, T-junction-free, and pure-quad internal tessellation of its patches. This requires the insertion of additional irregular-vertices inside patches, but the regularity of the final- mesh is safeguarded by optimizing for both their number and for their reciprocal alignment. In total, our method guarantees the reproduction of feature-lines by construction, while still producing good quality, isometric, pure-quad, conforming meshes, making it an ideal candidate for CAD models. Moreover, the method is fully automatic, requiring no user intervention, and remarkably reliable, requiring little assumptions on the input mesh, as we demonstrate by batch processing the entire Thingi10K repository, with less than 0.5% of the attempted cases failing to produce a usable mesh.

[Web Site](https://www.quadmesh.cloud/)
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
git clone --recursive https://github.com/TODO
```

### Compute Paper Results TODO

The code of this paper has been tested both on Linux Ubuntu 20.04 and MacOS.
To compute the results shown in the paper, you can run the following line on a shell:

```
sh scripts/[your platform]/make_it.sh
```

The script will download and install the required dependencies (cmake and cgal), it will build the code in a directory named `build`, and it will run the algorithm with the meshes contained in `misc/input_meshes` and with the proper parameters. Results (in the form of ply files) will be saved in a directory named `results`.

### Build TODO

To build the code, just run the following lines in a shell:

```
mkdir build
cd build
cmake ..
```

make sure to have installed cmake and CGAL in your machine. Inside the `build` directory a binary called `fourAxisMilling` will be created.

### Run TODO

The algorithm requires at least an input mesh file as parameter. Other parameters are not mandatory and default values are listed below.
To run the algorithm, from the `build` directory, just type in a shell:

```
./fourAxisMilling -i=input_mesh.obj [-o=output_directory] [parameters]
```

where parameters can be:

- `model_height`: the height of the output model along the rotation axes, expressed in millimeters; default value: 60;
- `stock_length`: the lenght of the raw cylinder stock, expressed in millimeters; default value: 100;
- `stock_diameter`: the diameter of the raw cylinder stock, expressed in millimeters; default value: 60;
- `prefiltering_smooth_iters`: number of Taubing smoothing iterations applied during prefiltering; default value: 500;
- `n_best_axis_dirs`: number of candidate direction for finding the best axis; default value: 2000;
- `n_visibility_dirs`: number of uniformly distributed visibility directions orthogonal of the rotation axis; default value: 120;
- `saliency_factor`: the saliency factor used for finding the segmentation using the graph-cut algorithm; default value: 25.0;
- `compactness_term`: the compactness term used for finding the segmentation using the graph-cut algorithm; default value: 30.0;
- `wall_angle`: angle between walls and fabrication direction; default value: 25.0;
- `max_first`: if this parameter is present, the block with +X direction will be considered as first block between top and bottom regions;
- `dont_scale_model`: if this parameter is present; the input mesh will not be scaled to fit into the stock;
- `just_segmentation`: if this parameter is present, the fabrication sequence (and the stocks-result shapes) won't be computed.

Some examples of runs:

```
./fourAxisMilling -i=kitten.obj -o=kitten_res --model_height=70 --stock_length=88 --stock_diameter=72

./fourAxisMilling -i=buddha.obj -o=buddha_res --model_height=70 --stock_diameter=72 --stock_length=86 --prefiltering_smooth_iters=750
```

## License
[GPL3](LICENSE) licensed
([FAQ](https://www.gnu.org/licenses/gpl-faq.html))



