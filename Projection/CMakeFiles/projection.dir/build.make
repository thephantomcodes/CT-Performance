# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dstring/code/CT-Performance/CT-Performance/Projection

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dstring/code/CT-Performance/CT-Performance/Projection

# Include any dependencies generated for this target.
include CMakeFiles/projection.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/projection.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/projection.dir/flags.make

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o: CMakeFiles/projection.dir/flags.make
CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o: src/ProjectionParameters.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dstring/code/CT-Performance/CT-Performance/Projection/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o -c /home/dstring/code/CT-Performance/CT-Performance/Projection/src/ProjectionParameters.cpp

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/projection.dir/src/ProjectionParameters.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dstring/code/CT-Performance/CT-Performance/Projection/src/ProjectionParameters.cpp > CMakeFiles/projection.dir/src/ProjectionParameters.cpp.i

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/projection.dir/src/ProjectionParameters.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dstring/code/CT-Performance/CT-Performance/Projection/src/ProjectionParameters.cpp -o CMakeFiles/projection.dir/src/ProjectionParameters.cpp.s

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.requires:

.PHONY : CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.requires

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.provides: CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.requires
	$(MAKE) -f CMakeFiles/projection.dir/build.make CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.provides.build
.PHONY : CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.provides

CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.provides.build: CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o


CMakeFiles/projection.dir/src/projection.cpp.o: CMakeFiles/projection.dir/flags.make
CMakeFiles/projection.dir/src/projection.cpp.o: src/projection.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dstring/code/CT-Performance/CT-Performance/Projection/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/projection.dir/src/projection.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/projection.dir/src/projection.cpp.o -c /home/dstring/code/CT-Performance/CT-Performance/Projection/src/projection.cpp

CMakeFiles/projection.dir/src/projection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/projection.dir/src/projection.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dstring/code/CT-Performance/CT-Performance/Projection/src/projection.cpp > CMakeFiles/projection.dir/src/projection.cpp.i

CMakeFiles/projection.dir/src/projection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/projection.dir/src/projection.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dstring/code/CT-Performance/CT-Performance/Projection/src/projection.cpp -o CMakeFiles/projection.dir/src/projection.cpp.s

CMakeFiles/projection.dir/src/projection.cpp.o.requires:

.PHONY : CMakeFiles/projection.dir/src/projection.cpp.o.requires

CMakeFiles/projection.dir/src/projection.cpp.o.provides: CMakeFiles/projection.dir/src/projection.cpp.o.requires
	$(MAKE) -f CMakeFiles/projection.dir/build.make CMakeFiles/projection.dir/src/projection.cpp.o.provides.build
.PHONY : CMakeFiles/projection.dir/src/projection.cpp.o.provides

CMakeFiles/projection.dir/src/projection.cpp.o.provides.build: CMakeFiles/projection.dir/src/projection.cpp.o


# Object files for target projection
projection_OBJECTS = \
"CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o" \
"CMakeFiles/projection.dir/src/projection.cpp.o"

# External object files for target projection
projection_EXTERNAL_OBJECTS =

projection: CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o
projection: CMakeFiles/projection.dir/src/projection.cpp.o
projection: CMakeFiles/projection.dir/build.make
projection: CMakeFiles/projection.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dstring/code/CT-Performance/CT-Performance/Projection/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable projection"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/projection.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/projection.dir/build: projection

.PHONY : CMakeFiles/projection.dir/build

CMakeFiles/projection.dir/requires: CMakeFiles/projection.dir/src/ProjectionParameters.cpp.o.requires
CMakeFiles/projection.dir/requires: CMakeFiles/projection.dir/src/projection.cpp.o.requires

.PHONY : CMakeFiles/projection.dir/requires

CMakeFiles/projection.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/projection.dir/cmake_clean.cmake
.PHONY : CMakeFiles/projection.dir/clean

CMakeFiles/projection.dir/depend:
	cd /home/dstring/code/CT-Performance/CT-Performance/Projection && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dstring/code/CT-Performance/CT-Performance/Projection /home/dstring/code/CT-Performance/CT-Performance/Projection /home/dstring/code/CT-Performance/CT-Performance/Projection /home/dstring/code/CT-Performance/CT-Performance/Projection /home/dstring/code/CT-Performance/CT-Performance/Projection/CMakeFiles/projection.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/projection.dir/depend

