# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.5.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.5.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jonathanhahn/Thesis/Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jonathanhahn/Thesis/Code/build

# Include any dependencies generated for this target.
include CMakeFiles/library.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/library.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/library.dir/flags.make

CMakeFiles/library.dir/src/core/rk45.c.o: CMakeFiles/library.dir/flags.make
CMakeFiles/library.dir/src/core/rk45.c.o: ../src/core/rk45.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/library.dir/src/core/rk45.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/library.dir/src/core/rk45.c.o   -c /Users/jonathanhahn/Thesis/Code/src/core/rk45.c

CMakeFiles/library.dir/src/core/rk45.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/library.dir/src/core/rk45.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/jonathanhahn/Thesis/Code/src/core/rk45.c > CMakeFiles/library.dir/src/core/rk45.c.i

CMakeFiles/library.dir/src/core/rk45.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/library.dir/src/core/rk45.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/jonathanhahn/Thesis/Code/src/core/rk45.c -o CMakeFiles/library.dir/src/core/rk45.c.s

CMakeFiles/library.dir/src/core/rk45.c.o.requires:

.PHONY : CMakeFiles/library.dir/src/core/rk45.c.o.requires

CMakeFiles/library.dir/src/core/rk45.c.o.provides: CMakeFiles/library.dir/src/core/rk45.c.o.requires
	$(MAKE) -f CMakeFiles/library.dir/build.make CMakeFiles/library.dir/src/core/rk45.c.o.provides.build
.PHONY : CMakeFiles/library.dir/src/core/rk45.c.o.provides

CMakeFiles/library.dir/src/core/rk45.c.o.provides.build: CMakeFiles/library.dir/src/core/rk45.c.o


CMakeFiles/library.dir/src/core/outputJson.c.o: CMakeFiles/library.dir/flags.make
CMakeFiles/library.dir/src/core/outputJson.c.o: ../src/core/outputJson.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/library.dir/src/core/outputJson.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/library.dir/src/core/outputJson.c.o   -c /Users/jonathanhahn/Thesis/Code/src/core/outputJson.c

CMakeFiles/library.dir/src/core/outputJson.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/library.dir/src/core/outputJson.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/jonathanhahn/Thesis/Code/src/core/outputJson.c > CMakeFiles/library.dir/src/core/outputJson.c.i

CMakeFiles/library.dir/src/core/outputJson.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/library.dir/src/core/outputJson.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/jonathanhahn/Thesis/Code/src/core/outputJson.c -o CMakeFiles/library.dir/src/core/outputJson.c.s

CMakeFiles/library.dir/src/core/outputJson.c.o.requires:

.PHONY : CMakeFiles/library.dir/src/core/outputJson.c.o.requires

CMakeFiles/library.dir/src/core/outputJson.c.o.provides: CMakeFiles/library.dir/src/core/outputJson.c.o.requires
	$(MAKE) -f CMakeFiles/library.dir/build.make CMakeFiles/library.dir/src/core/outputJson.c.o.provides.build
.PHONY : CMakeFiles/library.dir/src/core/outputJson.c.o.provides

CMakeFiles/library.dir/src/core/outputJson.c.o.provides.build: CMakeFiles/library.dir/src/core/outputJson.c.o


CMakeFiles/library.dir/src/core/parson/parson.c.o: CMakeFiles/library.dir/flags.make
CMakeFiles/library.dir/src/core/parson/parson.c.o: ../src/core/parson/parson.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/library.dir/src/core/parson/parson.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/library.dir/src/core/parson/parson.c.o   -c /Users/jonathanhahn/Thesis/Code/src/core/parson/parson.c

CMakeFiles/library.dir/src/core/parson/parson.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/library.dir/src/core/parson/parson.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/jonathanhahn/Thesis/Code/src/core/parson/parson.c > CMakeFiles/library.dir/src/core/parson/parson.c.i

CMakeFiles/library.dir/src/core/parson/parson.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/library.dir/src/core/parson/parson.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/jonathanhahn/Thesis/Code/src/core/parson/parson.c -o CMakeFiles/library.dir/src/core/parson/parson.c.s

CMakeFiles/library.dir/src/core/parson/parson.c.o.requires:

.PHONY : CMakeFiles/library.dir/src/core/parson/parson.c.o.requires

CMakeFiles/library.dir/src/core/parson/parson.c.o.provides: CMakeFiles/library.dir/src/core/parson/parson.c.o.requires
	$(MAKE) -f CMakeFiles/library.dir/build.make CMakeFiles/library.dir/src/core/parson/parson.c.o.provides.build
.PHONY : CMakeFiles/library.dir/src/core/parson/parson.c.o.provides

CMakeFiles/library.dir/src/core/parson/parson.c.o.provides.build: CMakeFiles/library.dir/src/core/parson/parson.c.o


CMakeFiles/library.dir/src/core/systems.c.o: CMakeFiles/library.dir/flags.make
CMakeFiles/library.dir/src/core/systems.c.o: ../src/core/systems.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/library.dir/src/core/systems.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/library.dir/src/core/systems.c.o   -c /Users/jonathanhahn/Thesis/Code/src/core/systems.c

CMakeFiles/library.dir/src/core/systems.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/library.dir/src/core/systems.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/jonathanhahn/Thesis/Code/src/core/systems.c > CMakeFiles/library.dir/src/core/systems.c.i

CMakeFiles/library.dir/src/core/systems.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/library.dir/src/core/systems.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/jonathanhahn/Thesis/Code/src/core/systems.c -o CMakeFiles/library.dir/src/core/systems.c.s

CMakeFiles/library.dir/src/core/systems.c.o.requires:

.PHONY : CMakeFiles/library.dir/src/core/systems.c.o.requires

CMakeFiles/library.dir/src/core/systems.c.o.provides: CMakeFiles/library.dir/src/core/systems.c.o.requires
	$(MAKE) -f CMakeFiles/library.dir/build.make CMakeFiles/library.dir/src/core/systems.c.o.provides.build
.PHONY : CMakeFiles/library.dir/src/core/systems.c.o.provides

CMakeFiles/library.dir/src/core/systems.c.o.provides.build: CMakeFiles/library.dir/src/core/systems.c.o


CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o: CMakeFiles/library.dir/flags.make
CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o: ../src/variationalPoincare/newtonsMethodforQuadraticSystem.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o   -c /Users/jonathanhahn/Thesis/Code/src/variationalPoincare/newtonsMethodforQuadraticSystem.c

CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/jonathanhahn/Thesis/Code/src/variationalPoincare/newtonsMethodforQuadraticSystem.c > CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.i

CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/jonathanhahn/Thesis/Code/src/variationalPoincare/newtonsMethodforQuadraticSystem.c -o CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.s

CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.requires:

.PHONY : CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.requires

CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.provides: CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.requires
	$(MAKE) -f CMakeFiles/library.dir/build.make CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.provides.build
.PHONY : CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.provides

CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.provides.build: CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o


# Object files for target library
library_OBJECTS = \
"CMakeFiles/library.dir/src/core/rk45.c.o" \
"CMakeFiles/library.dir/src/core/outputJson.c.o" \
"CMakeFiles/library.dir/src/core/parson/parson.c.o" \
"CMakeFiles/library.dir/src/core/systems.c.o" \
"CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o"

# External object files for target library
library_EXTERNAL_OBJECTS =

liblibrary.a: CMakeFiles/library.dir/src/core/rk45.c.o
liblibrary.a: CMakeFiles/library.dir/src/core/outputJson.c.o
liblibrary.a: CMakeFiles/library.dir/src/core/parson/parson.c.o
liblibrary.a: CMakeFiles/library.dir/src/core/systems.c.o
liblibrary.a: CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o
liblibrary.a: CMakeFiles/library.dir/build.make
liblibrary.a: CMakeFiles/library.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jonathanhahn/Thesis/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C static library liblibrary.a"
	$(CMAKE_COMMAND) -P CMakeFiles/library.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/library.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/library.dir/build: liblibrary.a

.PHONY : CMakeFiles/library.dir/build

CMakeFiles/library.dir/requires: CMakeFiles/library.dir/src/core/rk45.c.o.requires
CMakeFiles/library.dir/requires: CMakeFiles/library.dir/src/core/outputJson.c.o.requires
CMakeFiles/library.dir/requires: CMakeFiles/library.dir/src/core/parson/parson.c.o.requires
CMakeFiles/library.dir/requires: CMakeFiles/library.dir/src/core/systems.c.o.requires
CMakeFiles/library.dir/requires: CMakeFiles/library.dir/src/variationalPoincare/newtonsMethodforQuadraticSystem.c.o.requires

.PHONY : CMakeFiles/library.dir/requires

CMakeFiles/library.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/library.dir/cmake_clean.cmake
.PHONY : CMakeFiles/library.dir/clean

CMakeFiles/library.dir/depend:
	cd /Users/jonathanhahn/Thesis/Code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jonathanhahn/Thesis/Code /Users/jonathanhahn/Thesis/Code /Users/jonathanhahn/Thesis/Code/build /Users/jonathanhahn/Thesis/Code/build /Users/jonathanhahn/Thesis/Code/build/CMakeFiles/library.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/library.dir/depend

