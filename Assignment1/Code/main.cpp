#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

/**
 * @brief Computes the view matrix for a given eye position.
 *
 * This function computes the view matrix for a camera positioned at the
 * specified eye position. The view matrix is responsible for transforming
 * world coordinates to camera coordinates.
 *
 * @param eye_pos The position of the camera in world coordinates.
 * @return The view matrix representing the transformation from
 *         world to camera coordinates.
 */
Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0],
                 0, 1, 0, -eye_pos[1],
                 0, 0, 1, -eye_pos[2],
                 0, 0, 0, 1;

    view = translate * view;

    return view;
}

/**
 * @brief Computes the model matrix for a given rotation angle.
 *
 * This function computes the model matrix for rotating an object
 * around the Z axis by the specified rotation angle in degree.
 *
 * @param  rotation_angle The angle of rotation around the Z axis in radians.
 * @return The model matrix representing the rotation transformation.
 */
Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    // Convert the rotation angle from degree to radians
    float theta = rotation_angle * M_PI / 180.0f; 

    // Initialize the model matrix as identity matrix
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);

    // Modify the rotation part of the model matrix
    model.block<2, 2>(0, 0) << cos_theta, -sin_theta,
                               sin_theta, cos_theta;
    return model;
}

/**
 * @brief Computes the projection matrix for the given parameters.
 *
 * This function computes the projection matrix based on the given parameters
 * describing the view frustum.
 *
 * @param eye_fov The vertical field of view angle in degrees.
 * @param aspect_ratio The aspect ratio of the viewport (width/height).
 * @param zNear The distance to the near clipping plane.
 * @param zFar The distance to the far clipping plane.
 * @return The projection matrix representing the perspective projection.
 */
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear,   float zFar)
{
    // Variables to describe the vision field
    float l, r;    // Left, right
    float b, t;    // bottom, top
    float n, f;    // near, far

    // Compute the fov in radian and calculate half of the vertical fov angle
    float fov_rad = eye_fov * M_PI / 180.0f;
    float tan_half_fov = tan(fov_rad / 2.0f);

    // Calculate the param of cuboid
    t = zNear * tan_half_fov;
    b = -t;
    r = aspect_ratio * t;
    l = -r;
    n = zNear;
    f = zFar;

    // Calculate the projection matrix for perspective projection
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // Matrix for orthographic translation
    Eigen::Matrix4f ortho = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f ortho_scale = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f ortho_trans = Eigen::Matrix4f::Identity();
    ortho_trans << 1, 0, 0, -(r + l)/2.0f,
                   0, 1, 0, -(t + b)/2.0f,
                   0, 0, 1, -(n + f)/2.0f,
                   0, 0, 0, 1;
    ortho_scale << 2.0f / (r - l), 0,              0,              0,
                   0,              2.0f / (t - b), 0,              0,
                   0,              0,              2.0f / (n - f), 0,
                   0,              0,              0,              1;
    ortho = ortho_scale * ortho_trans;

    // Calculate the projection 2 orthographic matrix
    float A, B;    // A and B are the param in projection 2 ortho matrix
    A = n + f;
    B = -n * f;
    Eigen::Matrix4f proj2ortho = Eigen::Matrix4f::Identity();
    proj2ortho << n, 0, 0, 0,
                  0, n, 0, 0,
                  0, 0, A, B,
                  0, 0, 1, 0;

    // Now we calculate the projection matrix by multiple 2 matrix
    // M(orthographic) * M(projection 2 orthographic)
    projection = ortho * proj2ortho;

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
