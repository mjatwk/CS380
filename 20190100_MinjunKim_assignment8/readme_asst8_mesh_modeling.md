# About Mesh in Asst8
### Assignment 8 of 2023 Spring CS380: Introduction to Computer Graphics - 20190100 Minjun Kim

## Files
- model.mesh: identical to christmas_bear_2500.mesh
- mesh/
  - christmas_bear_{original|500|2500|5000}.{off|mesh}
  - head_{original|500|2500|5000}.{off|mesh}
  - off_to_mesh.py: python program to convert .off file into .mesh file
  - bear_2500.mov


## Converting and Decimation
Original files downloaded from [Artec3D](https://www.artec3d.com/3d-models) as .obj file. 

Used [MeshLab](https://www.meshlab.net/) to convert .obj file into .off file. 

Since the conditions were to be able to subdivide the mesh at least 3 levels, which made me to decimate the models. My machine (Macbook Air 2nd Gen Silicon Model) could only control up to ~300k faces, initial mesh could have maximum ~5k faces. But to enable the model to decrease the subdivision level back to 2, starting with 2.5k faces seemed to be rational for machine stability. In such reason, the 2.5k face mesh is selected as model.mesh, which versions with different number of faces are included in /mesh folder along with the original file as well.

Even though the .off file was similar enough with the .mesh file provided from the assignment, as the instructor has recommended, it did have some differences. This forced to make a converter from .off file to .mesh, which is named as **off_to_mesh.py**, also is submitted inside the /mesh folder.

## About the Submitted File
Please note that the bear submitted as model.mesh is lying on the z-coordinate, while watching the +x coordinate. For the face, it is facing the -y coordinate. Two videos head_2500.mov and bear_2500.mov is submitted to briefly show how the mesh works.