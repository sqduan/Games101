// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>

// Load the positions of vertices
rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

// Load the indices of vertices
rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

// Load the color of each vertices
rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

// Transfer from 
auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

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

#ifdef DEBUG
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "The triangle's vertices are: " << t.v[0].transpose() << " "
        << t.v[1].transpose() << " " << t.v[2].transpose() << std::endl;
    std::cout << "The bounding box vertices are: " << box.min_point.transpose() << " "
        << box.max_point.transpose() << std::endl;
#endif

    return; 
}

static bool insideTriangle(int x, int y, const Vector3f* _v)
{   
    bool is_inside = false;
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    Vector3f v0 = _v[1] - _v[0];
    Vector3f v1 = _v[2] - _v[1];
    Vector3f v2 = _v[0] - _v[2];

    float cross1 = (x - _v[0].x()) * v0.y() - (y - _v[0].y()) * v0.x();
    float cross2 = (x - _v[1].x()) * v1.y() - (y - _v[1].y()) * v1.x();
    float cross3 = (x - _v[2].x()) * v2.y() - (y - _v[2].y()) * v2.x();

    is_inside = (cross1 > 0 && cross2 > 0 && cross3 > 0) ||
                (cross1 < 0 && cross2 < 0 && cross3 < 0);
    return is_inside;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) /
        (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) /
        (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) /
        (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;

    // For each of the triangles
    for (auto& i : ind)
    {
        Triangle t;

        // Calculate the transformation for the vertices
        Eigen::Vector4f v[] = {
            mvp * to_vec4(buf[i[0]], 1.0f),
            mvp * to_vec4(buf[i[1]], 1.0f),
            mvp * to_vec4(buf[i[2]], 1.0f)
        };

        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }

        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        // Set the vertex of triangle on the screen
        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    bool is_inside;
    Eigen::Vector3f pixel;

    auto v = t.toVector4();
    
    // TODO : Find out the bounding box of current triangle.
    rst::bounding_box box;
    getTriangleBoundingBox(t, box);
    
    // iterate through the pixel and find if the current pixel is inside the triangle
    for (int y = int(box.min_point.y()); y < int(box.max_point.y()); y++) {
        for (int x = int(box.min_point.x()); x < int(box.max_point.x()); x++) {
            is_inside = insideTriangle(x, y, t.v);
            if (is_inside) {
                auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;

                if (z_interpolated < depth_buf[x * width + y]) {
                    pixel.x() = x;
                    pixel.y() = y;
                    pixel.z() = z_interpolated;

                    // Set the color of the pixel
                    set_pixel(pixel, t.getColor());
                    depth_buf[x * width + y] = z_interpolated;
                }

            }
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height - 1 - y) * width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on
