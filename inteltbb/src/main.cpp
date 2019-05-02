// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/resource.h>
#include <Eigen/Dense>
#include <time.h>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "shape.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

vector<Shape*> shape_vectors;
int num_thread;


void extend_stack() {
    const rlim_t kStackSize = 16 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    } else {
        fprintf(stderr, "getrlimit returned result = %d\n", result);
    }
}

int main(int argc, char *argv[])
{
    // clock_t tStart = clock();
    // The default stack limit is 8M; extend it to 16M just in case.
    extend_stack();

    if (argc != 3)
    {
        // std::cout << "Error: not enough input argument" << std::endl;
        printf("Error: not enough input argument");
        exit(1);
    }

    num_thread = atoi(argv[1]);

    string filename(argv[2]);



    Mesh *m = new Mesh(filename, LAMBERTIAN_SHADING, COLOR_GOLD);
    shape_vectors.push_back(m);



    compute_scene();

    for (vector<Shape*>::iterator it = shape_vectors.begin();
              it != shape_vectors.end(); ++it) {
      delete *it;
    }

    return 0;
}
