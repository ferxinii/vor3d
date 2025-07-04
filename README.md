# Voronoi 3D
This code allows the construction of bounded 3D Voronoi diagrams inside an arbitrary polytope. 

Given a set of seeds, the code first constructs the corresponding Delaunay triangulation using an iterative insertion flipping algorithm. Next, the dual Voronoi diagram is extracted. It is implemented in such a way that the Voronoi diagram is automatically bounded inside a specified polytope.


<p align="center">
<img src="./images/sph_v1.png" alt="Example sphere 1" width="400" height="auto" />
<img src="./images/sph_v3.png" alt="Example sphere 2" width="400" height="auto">
</p>
<p align="center">
<img src="./images/cube_v1.png" alt="Example cube 1" width="400" height="auto" />
<img src="./images/cube_v4.png" alt="Example cube 2" width="400" height="auto">
</p>


## Complex example
In the following example, we see how the code has been used to generate secondary lung lobules inside a lung geometry.

<p align="center">
<img src="./images/R_v2.png" alt="Example right lung" width="800" height="auto">
</p>
<p align="center">
<img src="./images/L_v3.png" alt="Example left lung" width="800" height="auto">
</p>
