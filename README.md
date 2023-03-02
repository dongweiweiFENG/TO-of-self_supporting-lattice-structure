# TO-of-self_supporting-lattice-structure
## Overview
![image](https://user-images.githubusercontent.com/124340386/222349655-d4113531-157c-4424-a1c8-a7f7aa017ab0.png)

![image](https://user-images.githubusercontent.com/124340386/222351163-927c871b-0e23-46e8-bafc-71dda1a7d818.png)
the density-based 3D topology optimization method (3DTop)  
the TO method with self-supporting filter (SS3DTop)  
the results of the ground lattice structure (GLS)  
the direct simplification on ground lattice structure(DSGLS)  

## Abstract
This paper presents a method to generate self-supporting lattice structures in a topology optimization (TO) framework.
Our method is composed of two phases, which are the TO-based subdivision and the TO-based simplification. Starting
from a lattice structure with self-supporting lattice unit cells in a coarse resolution, a subdivision method is proposed to
adaptively generate the lattice structure based on the density-based topology optimization framework. The subdivision
operators are well-designed to preserve self-supporting property on struts, and a filtering approach is proposed to avoid
overhanging nodes. To remove redundant struts on the lattice structures generated by subdivision, a simplification
method is developed to satisfy the required volume while a self-supporting constraint is incorporated to ensure the
manufacturability of the resultant structures. Our method has been tested on both 2D and 3D examples. Experimental
tests demonstrate that our method can effectively generate self-supporting lattice structures with stronger mechanical
strength.
## Usage
The code what we provide is for 2D examples which is easy to implement for ordinary Chixian personal notebook computer.
If you want to run the code, first make sure that MATLAB is installed on the computer. Then you need to download all .m files to your local folder. Through topOLAttice_ Subdivision .m file can realize the first step of segmentation and store the subdivision data in 2d_subdivision.mat file, and then use the topOLAttice_simplification .m file implements the second step of simplification.

## Contact Information
Weiming Wang(WWMdlut@gmail.com)  
Dongwei Feng(dutshuxuefeng@163.com)  
Charlie C.L. Wang (changling.wang@manchester.ac.uk)
