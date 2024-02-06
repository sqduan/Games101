# 作业2思路
在作业1中我们实现了对物体顶点的投射变换，在作业2中，我们要对三角形进行光栅化，从而绘制一个具有颜色的实心三角形，步骤如下：
1. 创建三角形的 2 维bounding box。
2. 遍历此 bounding box 内的所有像素（使用其整数索引）。然后，使用像素中心的屏幕空间坐标来检查中心点是否在三角形内。
3. 如果在内部，则将其位置处的插值深度值 (interpolated depth value) 与深度缓冲区 (depth buffer) 中的相应值进行比较。
4. 如果当前点更靠近相机，请设置像素颜色并更新深度缓冲区 (depth buffer)

# 作业2代码
## 创建三角形的Bounding box
我们首先定义一个结构体描述三角形的Bounding Box，该结构包含bounding box左下角和右上角两个点
```C++
// The struct to describe 2D bounding box
struct bounding_box {
    Eigen::Vector2f min_point;
    Eigen::Vector2f max_point;
};
```
bounding box的计算过程较为简单，我们只需要遍历三角形的三个顶点，分别确定`x_min, y_min, x_max, y_max`。此处定义一个函数`get_triangle_bounding_box`，该函数接收待计算bounding box的三角形`t`，并返回对应的bounding box。
```C++
static void getTriangleBoundingBox(const Triangle& t, rst::bounding_box& box)
{
    // Compute the minimum and maximum points of the bounding box
    box.min_point.x() = std::numeric_limits<float>::max();
    box.min_point.y() = std::numeric_limits<float>::max();
    box.max_point.x() = std::numeric_limits<float>::lowest();
    box.max_point.y() = std::numeric_limits<float>::lowest();

    // Update the bounding box value
    for (auto& v : t.v) {
        box.min_point.x() = std::min(box.min_point.x(), v.x());
        box.min_point.y() = std::min(box.min_point.y(), v.y());
        box.max_point.x() = std::max(box.max_point.x(), v.x());
        box.max_point.y() = std::max(box.max_point.y(), v.y());
    }

    return; 
}
```
## 遍历bounding box中所有的像素点，并判断是否在三角形内
获得bounding box后，下一步是遍历其中所有的像素点，判断每个像素点是否在三角形内部。设三角形的三个顶点构成三条向量，分别为：
$$
\mathbf{AB}, \ \mathbf{BC},\ \mathbf{CA}
$$
设像素点的坐标为$P$，则$P$与三个顶点构成三条向量：
$$
\mathbf{AP}, \ \mathbf{BP},\ \mathbf{CP}
$$
若$\mathbf{AB}\times\mathbf{AP}$，$\mathbf{BC}\times\mathbf{BP}$，$\mathbf{CA}\times\mathbf{CP}$同向，该点在三角形内部。
由于屏幕坐标系是二维的，两个二维向量的叉乘


