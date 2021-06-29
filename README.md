# meshlab_manual
translation of menu bar and help of some sampling and reconstruction algorithm

# Content
* [Filters](#Filters)  
    * [Remeshing, Simplification and Reconstruction](#Remeshing)
        * [Simplification: Quadric Edge Collapse Decimation](#二次边折叠抽取)
        * [Surface Reconstruction: Ball Pivoting](#表面重建之球旋转)
    * [Normals, Curvatures and Orientation](#Normals)
    * [Sampling](#Sampling)
        * [Clustered Vertex Sampling](#聚类顶点采样)

        * [Disk Vertex Coloring](#磁盘顶点着色)

        * [Distance from Reference Mesh](#与参考网格的距离)

        * [Dust Accumulation](#灰尘堆积)

        * [Hausdorff Distance](#豪斯多夫距离)

        * [Mesh Element Sampling](#网格元素采样)

        * [Montecarlo Sampling](#蒙特卡洛采样)

        * [Point Cloud Simplification](#点云简化)

        * [Poisson-disk Sampling](#泊松盘采样)

        * [Regular Recursive Sampling](#常规递归采样)

        * [Stratified Triangle Sampling](#分层三角采样)

        * [Texel Sampling](#纹素采样)

        * [Vertex Attribute Transfer](#顶点属性转移)

        * [Volumetric Sampling](#体积采样)

        * [Voronoi SamplingVoronoi](#采样)

        * [Voronoi ScaffoldingVoronoi](#脚手架)

        * [Voronoi Vertex ColoringVoronoi](#顶点着色)
    
===============================================

# Filters

## Remeshing, Simplification and Reconstruction

### Simplification: Quadric Edge Collapse Decimation
### 二次边折叠抽取
*Simplify a mesh using a Quadric based Edge Collapse Strategy; better than clustering but slower*  
*使用基于二次曲面的边缘折叠策略简化网格； 比聚类好但更慢*

* **Tarfet number of faces**  
The desired final number of faces.
期望的最终的面片数

* **Percentage reduction(0..1)**  
If not zero, this parameter specifies the desired final size of the mesh as a percentage of the initial size.  
如果不为零，则此参数将所需的网格最终尺寸指定为初始尺寸的百分比。

* **Quality threshold**  
Quality threshold for penalizing bad shaped faces.
The value is in the range [0..1] 0 accept and kind of face(no penalties), 0.5 penalize faces with quality  
惩罚形状不佳的面片的质量阈值。  
该数值范围在[0,1]，0表示接受任意形状的面片（即无惩罚），0.5表示根据面片质量惩罚。

* **Preserve Boundary of the mesh**  
The simplification process tries to do not affect mesh boundaries during simplification  
简化程序尽量不影响简化过程中的网格边界

    * **Boundary Preserving Weight**
    The importance of the boundary during simplification. Default (1..0) means that the boundary has the same importance of the rest. Values greater than 1.0 raise boundary importance and has the effect of removing less vertices on the border.  
    简化过程中边界的重要性。 默认 (1..0) 表示边界与其余边界具有相同的重要性。 大于 1.0 的值会提高边界的重要性，并具有移除边界上更少顶点的效果。

* **Preserve Normal**  
Try to avoid face flipping effects and try to preserve the original orientation of the surface  
尽量避免面翻转效果，尽量保持曲面的原始方向

* **Preserve Topology**  
Avoid all the collapses that should cause a topology change in the mesh (like closing holes, squeezing handles, etc). If checked the genus of the mesh should stay unchanged.  
避免所有会导致网格拓扑变化的塌陷（如关闭孔、挤压手柄等）。 如果选中，网格的类应保持不变。

* **Optimal position of simplified vertices**  
Each collapsed vertex is placed in the position minimizing the quadric error. It can fail (creating bad spikes) in case of very flat areas. If disabled edges are collapsed onto one of the two original vertices and the final mesh is composed.
每个折叠的顶点都放置在最小化二次误差的位置。 如果区域非常平坦，它可能会失败（产生不良尖峰）。 如果禁用，边折叠到两个原始顶点之一上，则最终网格就构成了。

* **Planar Simplification**  
Add additional simplification constraints that improves the quality of the simplification of the planar portion of the mesh, as a side effect, more triangles will be preseved in flat areas (allowing better shaped triangles).
添加额外的简化约束以提高网格平面部分的简化质量，作为副作用，将在平坦区域中保留更多三角形（允许形状更好的三角形）。

    * **Planar Simp. Weight**  
    How much we should try to preserve the triangles in the planar regions. If you lower this value planar areas will be simplified more.  
    我们应该尝试在平面区域中保留多少三角形。 如果您降低此值，平面区域将更加简化。

* **Weighted Simplification**  
Use the Per-Vertex quality as a weighting factor for the simplication value, so a vertex with a high quality value will not be simplified and a portion of the mesh with low quality values will be  
使用 Per-Vertex 质量作为简化值的加权因子，因此质量值高的顶点不会被简化，而质量值低的网格部分将被简化

* **Post-simplification cleaning**  
After the simplification an additional set of steps is performed to clean the mesh (unreferenced vertices, bad faces, etc)  
简化后，执行一组额外的步骤来清理网格（未引用的顶点、坏面等）

* **Simplify only selected faces**  
The simplification is applied only to the selected set of faces. Take care of the target number of faces!
简化仅应用于选定的一组面。 照顾目标数量的面孔！

### Surface Reconstruction: Ball Pivoting
### 表面重建之球旋转
*Given a point cloud with normals it reconstructs a surface using the Ball Pivoting Algorithm. Starting with a seed triangle, the BPA algorithm pivots a ball of the given radius around the already formed edges until it touches another point, forming another triangle. The process continues until all reachable edges have been tried. This surface reconstruction algorithm uses the existing points without creating new ones. Works better with uniformly sampled point cloud. I needed first perform a poisson disk subsampling of the point cloud.*  
*Bernardini F., Mittleman J., Rushmeier H., Silva C., Taubin G., **The ball-pivoting algorithm for surface reconstruction.***  
*IEEE TVCG 1999*  
*给定一个带有法线的点云，它使用球旋转算法重建一个表面。 从种子三角形开始，BPA 算法围绕已经形成的边缘旋转给定半径的球，直到它接触另一个点，形成另一个三角形。 这个过程一直持续到所有可到达的边都被尝试过。 此表面重建算法使用现有点而不创建新点。 使用均匀采样的点云效果更好。 我需要首先对点云执行泊松盘二次采样。*  

* Pivoting Ball radius(0 autoguess) (abs and %)   
*world unit & perc on(0 .. 0.0573658)*  
The radius of the ball pivoting (rolling) over the set of points. Gaps that are large than the ball radius will not be filled; similarly the smasll pits that are smaller than the ball radius will be filled.  
球在一组点上旋转（滚动）的半径。 大于球半径的间隙将不会被填充； 同样，小于球半径的小坑将被填充。  

* Clustering radius (% of ball radius)  
To avoid the creation of too small triangles, if a vertex is found too close to a previous one, it is clustered/merged with it.  
为避免创建太小的三角形，如果发现顶点与前一个顶点太近，则将其聚类/合并。 

* Angle Threshold (degrees)  
If we encounter a crease angle that is too large we should stop the ball rolling  
如果我们遇到过大的折痕角度，我们应该停止滚球

* Delete initial set of faces  
if true all the initial faces of the mesh are deleted and the whole surface is rebuilt from scratch. Otherwise the current faces are used as a starting point. Useful if you run the algorithm multiple times with an increasing ball radius.  
如果为 true，则删除网格的所有初始面，并从头开始重建整个表面。 否则，当前面将用作起点。 如果您以增加的球半径多次运行该算法，则很有用。

## Normals, Curvatures and Orientation

### Compute curvature principal directions  
### 计算曲率主方向  

### Compute normals for point sets
### 计算点集的法线  

### Cross Field Creation
### 区域交叉创建

### Cut mesh along crease edges  
### 沿折痕边缘切割网格

### Discrete Curvatures
### 离散曲率

### Invert Faces Orientation
### 反转面方向  

### Matrix: Freeze Current Matrix
### 矩阵：冻结当前矩阵

### Matrix: Invert Current Matrix
### 矩阵：反转当前矩阵

### Matrix: Reset Current Matrix
### 矩阵：重置当前矩阵

### Matrix: Set from translation/rotation/scale
### 矩阵：从平移/旋转/缩放设置

### Matrix: Set/Copy Transformation
### 矩阵：设置/复制转换 

### Normalize Face Normals
### 标准化面部法线

### Normalize Vertex Normals
### 标准化顶点法线

### Per Vertex Normal Function
### 每顶点法线函数

### Re-Compute Face Normals
### 重新计算面法线

### Re-Computer Per-Polygon Face Normals
### 重新计算每多边形面法线

### Re-Computer Vertex Normals
### 重新计算顶点法线

### Re-Orient all faces coherently
### 连贯地重新定位所有面

### Re-Orient vertex normals using cameras
### 使用相机重新定向顶点法线

### Smooths normals on a points sets
### 平滑点集上的法线

### Transform: Align to Principal Axis
### 变换：与主轴对齐

### Transform: Flip and/or swap axis  
### 变换：翻转和/或交换轴

### Transform: Rotate  
### 变换：旋转

### Transform: Rotate to Fit to a plane
### 变换：旋转以适应平面

### Transform: Scale, Normalize
### 变换：缩放、标准化

### Transform: Translate, Center, set Origin
### 变换：平移、居中、设置原点  

## Sampling

### Clustered Vertex Sampling
### 聚类顶点采样

### Disk Vertex Coloring
### 磁盘顶点着色

### Distance from Reference Mesh
### 与参考网格的距离

### Dust Accumulation
### 灰尘堆积

### Hausdorff Distance
### 豪斯多夫距离

### Mesh Element Sampling
### 网格元素采样

### Montecarlo Sampling
### 蒙特卡洛采样

### Point Cloud Simplification
### 点云简化

### Poisson-disk Sampling
### 泊松盘采样
*Create a new layer populated with a point sampling of the current mesh; samples are generated according to a Poisson-disk distribution, using the algorithm described in:*  
*创建一个新层，填充当前网格的点采样； 样本是根据泊松盘分布生成的，使用的算法如下：*  
***'Efficient and Flexible Sampling with Blue Noise Properties of Triangular Meshes'***  
*Massimiliano Corsini, Paolo Cignoni, Roberto Scopigno*  
*IEEE TVCG 2012*  

* Number of samples  
The desired number of samples. The ray of the disk is calculated according to the sampling density.  
所需的样本数。 根据采样密度计算圆盘的射线。  

* Explicit Radius (abs and %)  
world unit & perc on (0 .. 0.0571015)  
If not zero this parameter override the previous parameter to allow exact radius specification  
如果不为零，此参数将覆盖先前的参数以允许精确的半径规范

* MonterCarlo OverSampling  
The over-sampling rate that is used to generate the initial Montecarlo samples (e.g. if this parameter is K means that K x poisson sample points will be used). The generated Poisson-disk samples are a subset of these initial Montecalo samples. Larger this number slows the process but makes it a bit more accurate.  
用于生成初始蒙特卡洛样本的过采样率（例如，如果此参数为 K，则表示将使用 K x 泊松样本点）。 生成的泊松盘样本是这些初始 Montecalo 样本的子集。 更大的这个数字会减慢过程，但会使其更准确。

* Save Montecarlo  
If true, it will generate an additional Layer with the montecarlo sampling that was pruned to build the poisoon distribution.  
如果为真，它将生成一个额外的层，其中包含经过修剪以构建泊松分布的蒙特卡罗采样。  

* Approximate Geodesic Distance  
If true Poisson Disc distances are computed using an approximate geodesic distance, e.g. a euclidean distance weighted by a function of the difference between the normals of the two points.  
如果真正的泊松圆盘距离是使用近似测地距离计算的，例如 由两点法线之间的差值的函数加权的欧几里德距离。  

* Base Mesh Subsampling  
If true the original vertices of the base mesh are used as base set of points. In this case the SampleNum should be obviously much smaller than the original vertex number.  
Note that this option is very useful in the case you want to subsample a dense point cloud.  
如果为 true，则基础网格的原始顶点将用作基础点集。 在这种情况下，SampleNum 显然应该比原来的顶点数小很多。  
请注意，此选项在您想要对密集点云进行二次采样的情况下非常有用。

* Refine Existing Samples  
If true the vertices of the below mesh are used as starting vertices, and they will utterly refined by adding more and more points until possible.  
如果为 true，则下面网格的顶点将用作起始顶点，并且通过添加越来越多的点直到可能为止，它们将完全细化。  

* Samples to be refined  
Used only if the above option is checked.  
仅在选中上述选项时使用。

* Best Sample Heuristic  
If true it will use a simple heuristic for choosing the samples. At a small cost (it can slow a bit the process) it usually improve the maximality of the generated sampling.  
如果为 true，它将使用简单的启发式方法来选择样本。 以很小的代价（它可以稍微减慢过程），它通常可以提高生成采样的最大值。  

* Best Sample Pool Size  
Used only if the Best Sample Flag is true. It control the number of attempt that it makes to get the best sample. It is reasonable that it is smaller than the Montecarlo oversampling factor.  
仅当最佳样本标志为真时使用。 它控制为获得最佳样本而进行的尝试次数。 它小于 Montecarlo 过采样因子是合理的。

* Exact number of samples  
If requested it will try to do a dichotomic search for the best poisson disk radius that will generate the requested number of samples with a tolerance of the 0.5%. Obviously it takes much longer.  
如果需要，它将尝试对最佳泊松盘半径进行二分搜索，以生成所需数量的样本，容差为 0.5%。 显然，这需要更长的时间。

* Radius Variance  
The radius of the disk is allowed to vary between r and r*var. If this parameter is 1 the sampling is the same of the Poisson Disk Sampling  
圆盘的半径允许在 r 和 r*var 之间变化。 如果此参数为 1，则采样与泊松盘采样相同  

### Regular Recursive Sampling
### 常规递归采样

### Stratified Triangle Sampling
### 分层三角采样

### Texel Sampling
### 纹素采样

### Vertex Attribute Transfer
### 顶点属性转移

### Volumetric Sampling
### 体积采样

### Voronoi Sampling
### Voronoi采样

### Voronoi Scaffolding
### Voronoi脚手架

### Voronoi Vertex Coloring
### Voronoi顶点着色

