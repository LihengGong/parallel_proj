//
// Created by LihengGong & Hao Huang on 03/20/19.
//
// Shape related data structure and method definition
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include "shape.h"
#include "utils.h"


using namespace std;
using namespace Eigen;
extern vector<Shape*> shape_vectors;

const double BLIN_PHONE_EXPONENTIAL = 20;


double color_matrix[][3] = {
  {1., 0.,     0.}, //RED
  {0., 1.,     0.}, //GREEN
  {0., 0.,     1.}, //BLUE
  {1., 0.8431, 0.}, //GOLD
  {196/255, 51/255, 1.}, //NYU
  {1., 1.,     1.}, //GREY
};

// Camera scene related data members.
// These members are common to all shapes.
const unsigned PIXELS = 600;
long Shape::x_len = PIXELS;
long Shape::y_len = PIXELS;
MatrixXd Shape::C_R_mat = MatrixXd::Zero(Shape::x_len, Shape::y_len);
MatrixXd Shape::C_G_mat = MatrixXd::Zero(Shape::x_len, Shape::y_len);
MatrixXd Shape::C_B_mat = MatrixXd::Zero(Shape::x_len, Shape::y_len);
MatrixXd Shape::BK_mat = MatrixXd::Ones(Shape::x_len, Shape::y_len);
MatrixXd Shape::A_mat = MatrixXd::Ones(Shape::x_len, Shape::y_len);

Vector3d Shape::pixel_origin(-10.0, 10.0, 4.0);
double Shape::scene_width = 20.0;
double Shape::scene_height = 20.0;
double Shape::camera_z_axis = 8.0;
bool** Shape::pixel_is_shading_bitmap = new bool*[Shape::x_len];
Vector3d Shape::scene_ray_origin(4.0, 4.0, 8.0);
Vector3d Shape::scene_ray_direction(0., 0., -1.0);
Vector3d Shape::orthview_direction(0., 0., -1.0);
Vector3d Shape::persview_origin(0., 0., 6.0);
Vector3d Shape::light_position(0., -4.0, 4.0);
Vector3d Shape::second_light_position(1.0, 1.0, 1.0);
int Shape::light_number = 1;
View Shape::view = ORTHOGRAPHIC_VIEW;
bool Shape::is_draw_shadow = false;
double Shape::background_color = 0.15;
Vector3d Shape::x_displacement(Shape::scene_width / Shape::x_len, 0., 0.);
Vector3d Shape::y_displacement(0., -Shape::scene_height / Shape::y_len, 0.);
Vector3d Shape::surface_normal(0., 0., 0.);
Vector3d Shape::final_intersection(0., 0., 0.);



const int MAX_RECURSIVE_DEPTH = 2;
Configuration compute_raycolor(const Ref<const Vector3d> point_vect,
         const Ref<const Vector3d> direction, int depth) {

  Configuration config;
  //////////////// reflection recursive call
  if (depth > MAX_RECURSIVE_DEPTH) {
    return config;
  }

  double nearest_t = std::numeric_limits<double>::infinity();
  bool is_intersect = false;
  Vector3d unit_normal(1., 0., 0.);
  Shading cur_shading;
  Color cur_color_enum;

  for (vector<Shape*>::iterator it = shape_vectors.begin();
            it != shape_vectors.end(); ++it) {

    double cur_solution_t;
    if ((*it)->hit_slow(point_vect, direction, cur_solution_t, unit_normal)) {
      if (cur_solution_t < nearest_t) {
        nearest_t = cur_solution_t;
        is_intersect = true;
        cur_shading = (*it)->shading;
        cur_color_enum = (*it)->color_enum;
      }
    }
  }

  if (is_intersect) {
    double cur_solution_t;
    //// intersection point
    //// light direction
    Vector3d intersection = point_vect + nearest_t * direction;
    //// shadow direction points from intersection to light source
    Vector3d shadow_direction = Shape::light_position - intersection;
#ifdef CG_SHADOW
    Vector3d start_point = intersection + SHADOW_EPSILON * shadow_direction;
    for (vector<Shape*>::iterator it = shape_vectors.begin();
              it != shape_vectors.end(); ++it) {
      bool is_shadow = (*it)->hit(start_point, shadow_direction, cur_solution_t);
      if (is_shadow) {
        return config;
      }
    }
#endif
    ////////// color
    double color_val = unit_normal.dot(shadow_direction.normalized());
    color_val = std::max(0., color_val);
    double pixel_value_new = 0.;

    if(cur_shading == BLINN_PHONG_SHADING) {
      Eigen::Vector3d vector_l = shadow_direction.normalized();
      Eigen::Vector3d vector_v = (point_vect - intersection).normalized();
      Eigen::Vector3d vector_h = (vector_l + vector_v).normalized();
      pixel_value_new = unit_normal.dot(vector_h);
      pixel_value_new = std::max(double(0.), pixel_value_new);
      pixel_value_new = std::pow(pixel_value_new, double(BLIN_PHONE_EXPONENTIAL));
    }

    config.color_enum = cur_color_enum;
    config.shading = cur_shading;
    config.col_value = color_val + pixel_value_new;
    return config;
  }
  return config;
}

const int PIC_NUM = 1;

void compute_scene() {
  double displacement = Shape::scene_width / PIC_NUM;

  // hardcore openmp
  Shape::light_position << -4.0 + displacement, 10.0, 4.0;

  #pragma omp parallel for schedule(dynamic)
  for (unsigned i = 0; i < Shape::x_len; i++)
  {
    /*** omp id ****/
    int threadid = omp_get_thread_num();
    if(threadid == 0 && i == 0){
      cout << "num of thread " << omp_get_num_threads() << endl;
    }

    for (unsigned j = 0; j < Shape::y_len; j++)
    {
      bool is_shading = false;
      Eigen::Vector3d scene_ray_origin = Shape::pixel_origin +
              double(i) * Shape::x_displacement +
              double(j) * Shape::y_displacement;

      Configuration config = compute_raycolor(scene_ray_origin,
                               Shape::orthview_direction, 0);

      ///// This code is ugly. I hate it.
      switch (config.color_enum)
      {
      case COLOR_RED:
        Shape::C_R_mat(i, j) = config.col_value;
        break;
      case COLOR_GREEN:
        Shape::C_G_mat(i, j) = config.col_value;
        break;
      case COLOR_BLUE:
        Shape::C_B_mat(i, j) = config.col_value;
        break;
      case COLOR_GOLD:
        Shape::C_R_mat(i, j) = config.col_value;
        Shape::C_G_mat(i, j) = config.col_value * 0.8431;
        //Shape::C_B_mat(i, j) = config.col_value;
        break;
      case COLOR_NYU:
        Shape::C_R_mat(i, j) = config.col_value * (87.0 / 255.0);
        Shape::C_G_mat(i, j) = config.col_value * (6.0 / 255.0);
        Shape::C_B_mat(i, j) = config.col_value * (140.0 / 255.0);
        break;
      case COLOR_GREY:
        Shape::C_R_mat(i, j) = config.col_value;
        Shape::C_G_mat(i, j) = config.col_value;
        Shape::C_B_mat(i, j) = config.col_value;
        break;
      default:
        break;
      }
    }
  }

  // Write pixel values into png file.
  std::string filename = "filename" + std::to_string(0) + ".png";
  write_matrix_to_png(Shape::C_R_mat, Shape::C_G_mat,
                        Shape::C_B_mat, Shape::A_mat,
                        filename);

}


void Mesh::load_mesh_data(const std::string& filename) {
    std::string line;
    std::ifstream infile;
    infile.open(filename);

    const int VERTICES_EACH_SHAPE = 3;

    if (! infile.eof()) {
        std::getline(infile, line);
        if (line != "OFF")
        {
            fprintf(stderr, "wrong file format.\n");
            return;
        }
    }

    if (! infile.eof()) {
        std::getline(infile, line);
        stringstream stream(line);
        if (! stream.eof()) {
            stream >> vert_number;
        } else {
            return;
        }
        if (! stream.eof()) {
            stream >> tri_number;
        } else {
            return;
        }
        if (! stream.eof()) {
            int temp;
            stream >> temp;
        } else {
            return;
        }
    }

    vertices_matrix = MatrixXd::Zero(vert_number, VERTICES_EACH_SHAPE);
    triangles_matrix = MatrixXi::Zero(tri_number, VERTICES_EACH_SHAPE);

    for (int line_number = 0; line_number < vert_number; line_number++) {
        getline(infile, line);
        stringstream stream(line);
        for (int index = 0; index < VERTICES_EACH_SHAPE; index++) {
            stream >> vertices_matrix(line_number, index);
        }
    }

    for (int line_number = 0; line_number < tri_number; line_number++) {
        getline(infile, line);
        stringstream stream(line);
        int temp;
        stream >> temp;
        for (int index = 0; index < VERTICES_EACH_SHAPE; index++) {
            stream >> triangles_matrix(line_number, index);
        }
    }

    if (filename == "../data/bunny.off") {
      cout << "scale bunny" << endl;
      vertices_matrix = 38 * vertices_matrix;
      vertices_matrix.col(0) += 5.8 * MatrixXd::Ones(vert_number, 1);
      vertices_matrix.col(1) += 2.2 * MatrixXd::Ones(vert_number, 1);
      // vertices_matrix.col(2) -= 2.0 * MatrixXd::Ones(vert_number, 1);
    } else {
      cout << "scale bumpy cube" << endl;
      vertices_matrix = 0.6 * vertices_matrix;
    }

    vertices_matrix.col(1) += 2.8 * MatrixXd::Ones(vert_number, 1);
}
