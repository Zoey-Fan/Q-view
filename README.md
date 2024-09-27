
## Q-view: An Efficient and Fine-Grained Viewshed Analysis method in a Three-Dimensional Urban Complex Environment

### 1. Introduction

Q-view method involves creating an indexing model for large-scale datasets and enabling parallel querying between the line-of-sight (LOS) and the model. The proposed Q-View method enables efficient
and spatially-exhaustive analysis, effectively mitigating the complexities associated with traditional viewshed computations. The core workflow of Q-view is mainly simplified to two parts: index model construction and 
viewshed query.

##### ●  Index Model Construction.

Given the immense number
Given the number of meshes in these models, which can reach billions within an area of 10 square kilometers, performing viewshed 
calculations on such a vast scale incurs significant computational overhead, resulting in low efficiency. To address this
issue, it is essential to pre-organize and manage data of such a large scale to minimize frequent access to irrelevant
data. We utilize 3D R-tree for modeling, enabling efficient querying. 

##### ● Viewshed Query:
Further transform the visibility analysis issue into a spatial query problem. 
Based on the index model constructed in the first step, perform parallel line-of-sight 
queries to check for the presence of obstructions, thereby determining visibility.


### 2. Settings
   
#### 2.1. Program Settings

CPU： 32 cores,Intel(R) Xeon(R) E5-2620@ 2.10GHz

memory: 512G

Operating System: Ubuntu 20.04
 
#### 2.2. Program Dependencies

* PROJ (recommended version 9.1.1)

```javascript
sudo apt-get install proj-bin libproj-dev
```

* CGAL (recommended version 5.5.2)

```javascript
cd '/home/Downloads/CGAL-5.5.2'    
./configure 
make			
make install
```

* mpirun (recommended version 3.3.2)

```javascript
sudo apt-get install mpich
```

### 3. Operating Steps

#### 3.1. Index Model Construction Phase

The `build3dindex` folder is used to build indexes for each `OBJ file`. In this context, test_dataset refers to the test dataset. Please download the four data files `tile_23_26_OBJ`, `tile_23_27_OBJ`, `tile_24_26_OBJ`, and `tile_24_27_OBJ` from the website https://portal.csdi.gov.hk/. After downloading, extract them and save them in the `build3dindex/test_dataset/paper_area500/` folder. The test_indexes folder stores the generated index files. There are two header files: `RTree.h` and `query3dRtree.h`. `RTree.h` file is used to generate a 3D R-tree index, based on code written by Greg Douglas (http://www.auran.com), to which we added our 3D viewshed query functionality. The implementation details are in the `query3dRtree.h` file. The main construction code is in `build3dindex.cpp`. Execute the following code to generate an index model for this batch of data, and store the results in the `test_indexes` file.

```javascript
cd build3dindex
make clean
make
conda activate
python3 run.py
```
#### 3.2. Index Model Merging Phase
In this stage, the index models generated in the first step are merged to obtain a 
complete index model for the area. 
Navigate to the `mergeindex` folder and execute the `Mergeindex` executable file. 
The input data should be the folder containing multiple index models, 
and the output data will be the index model for the entire area.

```javascript
cd mergeindex
make clean
make
./Mergeindex --inputdir ../build3dindex/test_indexes --output ./indexes/test_area.3idx
```

#### 3.3. Viewshed Query:
The `Q-view` folder contains all the code for viewshed analysis. It includes `Q-view.cpp` and `Q-curve.cpp`, which represent the operations for viewshed querying and line-of-sight querying along curves, respectively. 
The input data should include the observation range in meters, the observer's longitude and latitude 
coordinates separated by spaces, and the observation height. The output data will consist of `PNG` format images representing the square area centered on the observer with the 
observation range as the radius and a `CSV` file storing the visibility information.


```javascript
cd Q-view
make clean
make
./Q-VIEW --index ../mergeindex/indexes/test_area.3idx
```
For the code related to curve visibility analysis, it is stored in `Q-curve.cpp`. 
Execute the `Q-CURVE` executable file. The input data should include curve parameters, 
which are the start and end point coordinates in latitude and longitude, speed, and pitch angle 
information. The output data will be `CSV` files containing the total trajectory information, 
visible trajectory segment information, and non-visible trajectory segment information. 
By executing the `image.py` file, you can obtain a cumulative viewshed map.


```javascript
cd Q-view
make clean
make
./Q-CURVE --index ../mergeindex/indexes/test_area.3idx --positionStart 114.166 22.2825 1 --positionEnd 114.166 22.280 1 --velocity 50 --angle 50
```

### 4. Contact

* Institution

Yifan Zhang @ National University of Defense Technology
 
Mengyu Ma @ National University of Defense Technology

* Email

  zhangyifan.000@nudt.edu.cn
  
  mamengyu10@nudt.edu.cn









