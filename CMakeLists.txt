cmake_minimum_required(VERSION 3.11)
project(SNARL_PROJECT_CXX)

# Enable FetchContent module
include(FetchContent)

# Fetch Eigen dependency (this part you already have)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0  # Replace with the version you want
)

FetchContent_MakeAvailable(Eigen)

# Fetch Boost dependency (Boost version 1.77.0 is just an example)
FetchContent_Declare(
  Boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-1.77.0  # Replace with your desired version
)

FetchContent_MakeAvailable(Boost)

# Add source files
set(SRCS 
    src/matrix.cpp
    src/arg_parser.cpp
    src/snarl_parser.cpp
    src/vcf_parser.cpp
    src/binary_analysis.cpp
    src/main.cpp
    src/quantitative_analysis.cpp
)

# Create an executable from the source files
add_executable(${PROJECT_NAME} ${SRCS})

# Link Eigen and Boost to the executable
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen Boost::Boost)

# Include directories for Boost and Eigen
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
