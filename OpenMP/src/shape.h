#ifndef CG_HW1_SHAPE_H
#define CG_HW1_SHAPE_H



#define CG_SHADOW

#include <vector>

#include </usr/local/include/omp.h>

enum Shading{
    LAMBERTIAN_SHADING,
    BLINN_PHONG_SHADING,
    AMBIENT_SHADING,
    MIRROR_SHADING
};

enum View{
    ORTHOGRAPHIC_VIEW,
    PERSPECTIVE_VIEW,
};

enum Color {
  COLOR_RED,
  COLOR_GREEN,
  COLOR_BLUE,
  COLOR_GOLD,
  COLOR_NYU,
  COLOR_GREY,
  COLOR_BLACK
};

struct Configuration {
  Color color_enum;
  Shading shading;
  double col_value;
  explicit Configuration(Color c=COLOR_BLACK,
                Shading s=LAMBERTIAN_SHADING, double v=0.) :
  color_enum(c), shading(s), col_value(v) {}
};

const double BLACK_COLOR = 0.;

const double EPSILON = 0.0000001;
const double SHADOW_EPSILON = 0.6;

/*
Moler-Trumbore is the fastest one.
This method is adapted from the code in this thesis:
https://cadxfem.org/inf/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
*/
inline bool MolerTrumbore(const Eigen::Ref<const Eigen::Vector3d> rayOrigin,
                          const Eigen::Ref<const Eigen::Vector3d> rayVector,
                          const Eigen::Ref<const Eigen::Vector3d> vert0,
                          const Eigen::Ref<const Eigen::Vector3d> vert1,
                          const Eigen::Ref<const Eigen::Vector3d> vert2,
                          Eigen::Ref<Eigen::Vector3d> intersection,
                          double& solution_t) {

  Eigen::Vector3d edge1 = vert1 - vert0;
  Eigen::Vector3d edge2 = vert2 - vert0;

  Eigen::Vector3d pvec = rayVector.cross(edge2);
  double det = edge1.dot(pvec);

  if (det > - EPSILON && det < EPSILON) {
      return false;
  }

  double inv_det = 1 / det;
  Eigen::Vector3d tvec = rayOrigin - vert0;
  float u = (tvec.dot(pvec)) * inv_det;
  if (u < 0. || u > 1.0) {
      return false;
  }

  Eigen::Vector3d qvec = tvec.cross(edge1);
  double v = (rayVector.dot(qvec)) * inv_det;
  if (v < 0. || u + v > 1.0) {
      return false;
  }

  solution_t = (edge2.dot(qvec)) * inv_det;

  intersection = rayOrigin + rayVector * solution_t;
  return true;
}

void compute_scene();


class Shape {
public:
  Shape(Shading s, Color clr=COLOR_BLACK) {
    shading = s;
    unit_normal << 0., 0., 0.;
    color_enum = clr;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  static Eigen::MatrixXd C_R_mat;
  static Eigen::MatrixXd C_G_mat;
  static Eigen::MatrixXd C_B_mat;
  static Eigen::MatrixXd A_mat;
  static Eigen::MatrixXd BK_mat;
  static long x_len, y_len;
  static bool** pixel_is_shading_bitmap;
  static Eigen::Vector3d pixel_origin;
  static Eigen::Vector3d x_displacement;
  static Eigen::Vector3d y_displacement;
  static Eigen::Vector3d scene_ray_origin;
  static Eigen::Vector3d scene_ray_direction;
  static Eigen::Vector3d orthview_direction;
  static Eigen::Vector3d persview_origin;
  static double scene_width, scene_height;
  static double closest_pixel_in_ray;
  static double camera_z_axis;
  static double background_color;
  static Eigen::Vector3d light_position;
  static Eigen::Vector3d second_light_position;
  Eigen::Vector3d unit_normal;
  static Eigen::Vector3d surface_normal;
  static Eigen::Vector3d final_intersection;
  static int light_number;
  static View view;
  static bool is_draw_shadow;
  Color color_enum;
  Shading shading;

  static inline void generate_camera_rays(int x, int y) {
    Shape::scene_ray_origin = Shape::pixel_origin + double(x) * Shape::x_displacement +
                                double(y) * Shape::y_displacement;
    if (Shape::view == ORTHOGRAPHIC_VIEW) {
      Shape::scene_ray_direction = Shape::orthview_direction;
    } else {
      Shape::scene_ray_direction = Shape::scene_ray_origin - Shape::persview_origin;
      Shape::scene_ray_origin = Shape::persview_origin;
    }
  }

  virtual inline bool hit_slow(const Eigen::Ref<const Eigen::Vector3d> point_vec,
           const Eigen::Ref<const Eigen::Vector3d> direction,
           double& solution_t,
           Eigen::Ref<Eigen::Vector3d> unit_normal) = 0;


  virtual inline bool hit(const Eigen::Ref<const Eigen::Vector3d> point_vec,
           const Eigen::Ref<const Eigen::Vector3d> direction,
           double& solution_t) = 0;

  virtual inline bool render_shape(int x, int y) {return true;}

  virtual ~Shape() {}
};

class Mesh: public Shape {
public:
  Eigen::Vector3d p_vec;
  Eigen::Vector3d dir;
  explicit Mesh(const std::string& filename,
                Shading s=LAMBERTIAN_SHADING, Color clr=COLOR_BLACK)
  : Shape(s, clr) {
    load_mesh_data(filename);
    if (tri_number > 0) {
      solution_array = new double[tri_number];
      is_intersection_arry = new bool[tri_number];
      for (int ind = 0; ind < tri_number; ++ind) {
        solution_array[ind] = std::numeric_limits<double>::infinity();
        is_intersection_arry[ind] = false;
      }
    } else {
      solution_array = NULL;
    }
    p_vec << 0., 0., 0.;
    dir << 0., 0., 0.;
  }

  virtual ~Mesh() override {
    delete[] solution_array;
    delete[] is_intersection_arry;
    solution_array = NULL;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  Eigen::MatrixXd vertices_matrix;
  Eigen::MatrixXi triangles_matrix;
  double *solution_array;
  bool *is_intersection_arry;
  long vert_number, tri_number;
  void load_mesh_data(const std::string& filename);

  inline bool hit_slow(const Eigen::Ref<const Eigen::Vector3d> point_vec,
           const Eigen::Ref<const Eigen::Vector3d> direction,
           double& solution_t,
           Eigen::Ref<Eigen::Vector3d> unit_norm) override {

    Eigen::Vector3d intersection(0., 0., 0.);
    bool is_intersect= false;
    bool is_shading = false;
    double temp_solution = std::numeric_limits<double>::infinity();
    int final_index = -1;

    // #pragma omp parallel for schedule(static)
    for (unsigned int ind = 0; ind < tri_number; ind++) {
      Eigen::Vector3i cur_triangle_vertices = triangles_matrix.row(ind);
      Eigen::Vector3d vertice0 = vertices_matrix.row(cur_triangle_vertices(0));
      Eigen::Vector3d vertice1 = vertices_matrix.row(cur_triangle_vertices(1));
      Eigen::Vector3d vertice2 = vertices_matrix.row(cur_triangle_vertices(2));

      double cur_solution;
      static int cnt_once = 0;
      int threadid = omp_get_thread_num();
      if(threadid == 0 && cnt_once == 0){
        std::cout << "num of thread " << omp_get_num_threads() << std::endl;
        cnt_once = 1;
      }

      is_intersect = MolerTrumbore(point_vec,
                                        direction,
                                        vertice0, vertice1, vertice2,
                                        intersection, cur_solution);

      if (is_intersect && cur_solution < temp_solution) {
        temp_solution = cur_solution;
        if (temp_solution < 0) {
          std::cout << "waning! negative solution" << std::endl;
        }
        final_index = ind;
        is_shading = true;
      }
    }

    if (is_shading) {
      Eigen::Vector3i cur_triangle_vertices = triangles_matrix.row(final_index);
      Eigen::Vector3d vertice0 = vertices_matrix.row(cur_triangle_vertices(0));
      Eigen::Vector3d vertice1 = vertices_matrix.row(cur_triangle_vertices(1));
      Eigen::Vector3d vertice2 = vertices_matrix.row(cur_triangle_vertices(2));
      unit_norm = ((vertice1-vertice0).cross(vertice2-vertice0));
      unit_norm = unit_norm.normalized();

      if(unit_norm.dot(direction) > 0) {
        unit_norm = - unit_norm;
      }
      solution_t = temp_solution;
      return true;
    }

    return false;
  }

  inline bool render_shape(int x, int y) override {
    return render_all_triangles(x, y);
  }

  inline bool hit(const Eigen::Ref<const Eigen::Vector3d> point_vec,
           const Eigen::Ref<const Eigen::Vector3d> direction,
           double& nearest_solution_t) override {

    Eigen::Vector3d intersection(0., 0., 0.);
    for (unsigned t = 0; t < tri_number; t++) {
      Eigen::Vector3i cur_triangle_vertices = triangles_matrix.row(t);
      Eigen::Vector3d vertice0 = vertices_matrix.row(cur_triangle_vertices(0));
      Eigen::Vector3d vertice1 = vertices_matrix.row(cur_triangle_vertices(1));
      Eigen::Vector3d vertice2 = vertices_matrix.row(cur_triangle_vertices(2));

      double solution_t;

      bool is_intersect = MolerTrumbore(point_vec,
                                         direction,
                                         vertice0, vertice1, vertice2,
                                         intersection, solution_t);

      if (is_intersect && solution_t > EPSILON) {
          return true;
      }
    }

    return false;
  }

  inline bool render_all_triangles(int x, int y) {
    return true;
  }
};

#endif
