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
CMAKE_SOURCE_DIR = /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build

# Include any dependencies generated for this target.
include CMakeFiles/smooth_opengl3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/smooth_opengl3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/smooth_opengl3.dir/flags.make

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o: CMakeFiles/smooth_opengl3.dir/flags.make
CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o: /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/smooth_opengl3/smooth_opengl3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o   -c /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/smooth_opengl3/smooth_opengl3.c

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/smooth_opengl3/smooth_opengl3.c > CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.i

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/smooth_opengl3/smooth_opengl3.c -o CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.s

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.requires:

.PHONY : CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.requires

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.provides: CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.requires
	$(MAKE) -f CMakeFiles/smooth_opengl3.dir/build.make CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.provides.build
.PHONY : CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.provides

CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.provides.build: CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o


# Object files for target smooth_opengl3
smooth_opengl3_OBJECTS = \
"CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o"

# External object files for target smooth_opengl3
smooth_opengl3_EXTERNAL_OBJECTS =

bin/smooth_opengl3: CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o
bin/smooth_opengl3: CMakeFiles/smooth_opengl3.dir/build.make
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libGL.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libSM.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libICE.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libX11.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXext.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXrandr.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXi.so
bin/smooth_opengl3: lib/libglut.so.3.10.0
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libGL.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libSM.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libICE.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libX11.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXext.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXrandr.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
bin/smooth_opengl3: /usr/lib/x86_64-linux-gnu/libXi.so
bin/smooth_opengl3: CMakeFiles/smooth_opengl3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bin/smooth_opengl3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/smooth_opengl3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/smooth_opengl3.dir/build: bin/smooth_opengl3

.PHONY : CMakeFiles/smooth_opengl3.dir/build

CMakeFiles/smooth_opengl3.dir/requires: CMakeFiles/smooth_opengl3.dir/progs/demos/smooth_opengl3/smooth_opengl3.c.o.requires

.PHONY : CMakeFiles/smooth_opengl3.dir/requires

CMakeFiles/smooth_opengl3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/smooth_opengl3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/smooth_opengl3.dir/clean

CMakeFiles/smooth_opengl3.dir/depend:
	cd /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0 /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0 /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build /media/david/Hard_Disk/ENSTA/ANNEE_2/IN203/TP_self/Projet/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles/smooth_opengl3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/smooth_opengl3.dir/depend

