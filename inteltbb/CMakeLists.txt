cmake_minimum_required(VERSION 2.8.12)
project(parallel_proj)

### Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

if(APPLE)
set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/8.0.0/lib")
set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/8.0.0/include")
endif()

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(TBB_LINK_FLAGS "-ltbb")
  ###set(TBBROOT "/usr/local/Cellar/tbb/2018_U5")
  if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        message('match clang 1.')
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      message('match clang 2.' )
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      message('match clang 2.' $OpenMP_CXX)
    endif()
  ### set(TBBLIB_DEBUG "")
endif()

### Add src to the include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
LINK_DIRECTORIES(${LINK_DIRECTORIES} /opt/intel/compilers_and_libraries_2019.0.117/mac/tbb/lib)

find_library(TBB_LIB tbb)
if (TBB_LIB)
  message('Found_Intel_tbb_library.')
  set(TBB_LINK_FLAGS "-ltbb")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${TBB_LINK_FLAGS}")
endif()

include_directories("/usr/local/include" "/usr/local/opt/llvm/include")


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (OPENMP_FOUND)
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

### Include Eigen for linear algebra
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen")
### include_directories("${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb")
### Compile all the cpp files in src
file(GLOB SOURCES
"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(${PROJECT_NAME}_bin ${SOURCES})
link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
target_link_libraries(${PROJECT_NAME}_bin "${TBB_LIB}")

