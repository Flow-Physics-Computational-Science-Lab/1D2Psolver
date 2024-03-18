# For now I don't know the minumum version of cmake that I will be requiring.
# Otherwise, I will just need to check the mininum required for the other depen-
# dencies.
cmake_minimum_required(VERSION 3.10)

# Create a project
project(1D2Psolver VERSION 1.0)

# Make cmake check the C++ version required.
# Set the variable CMAKE_CXX_STANDARD to 11
#     the variable CMAKE_CXX_STANDARD_REQUIRED to True
#     the variable to output compile_files.json for nvim LSP
# Add debug flag.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Configure file?
#configure_file(TutorialConfig.h.in TutorialConfig.h)
#set(HEADER_FILES ../src/readers.h)

#add_library(readers readers.cpp ${HEADER_FILES})
#target_include_directories(readers PRIVATE 
#                           ../src 
#                           )

# Add files
file(GLOB sources ../src/*.cpp ../src/*.h)
file(GLOB schemes ../src/schemes/*.cpp ../src/schemes/*.h)

# Add an executable called Solver to the project
add_executable(1D2Psolver.x 
    ${sources} 
    ${schemes}
    )

# just for example add some compiler flags
#target_compile_options(1D2Psolver PUBLIC -Wall -Wfloat-conversion -g)

# I don't know how to use the target_include_directories to include ${PROJECT_BINARY_DIR}
#target_include_directories(2d_Solver PRIVATE 
#                           ../build 
#                           )
