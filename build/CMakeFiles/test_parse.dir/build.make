# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/fastF

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fastF/build

# Include any dependencies generated for this target.
include CMakeFiles/test_parse.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_parse.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_parse.dir/flags.make

CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o: CMakeFiles/test_parse.dir/flags.make
CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o: ../src/v1.0.9/fastq_parse.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fastF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o   -c /home/fastF/src/v1.0.9/fastq_parse.c

CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/fastF/src/v1.0.9/fastq_parse.c > CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.i

CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/fastF/src/v1.0.9/fastq_parse.c -o CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.s

CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o: CMakeFiles/test_parse.dir/flags.make
CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o: ../src/v1.0.9/argparse.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fastF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o   -c /home/fastF/src/v1.0.9/argparse.c

CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/fastF/src/v1.0.9/argparse.c > CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.i

CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/fastF/src/v1.0.9/argparse.c -o CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.s

# Object files for target test_parse
test_parse_OBJECTS = \
"CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o" \
"CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o"

# External object files for target test_parse
test_parse_EXTERNAL_OBJECTS =

test_parse: CMakeFiles/test_parse.dir/src/v1.0.9/fastq_parse.c.o
test_parse: CMakeFiles/test_parse.dir/src/v1.0.9/argparse.c.o
test_parse: CMakeFiles/test_parse.dir/build.make
test_parse: CMakeFiles/test_parse.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/fastF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable test_parse"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_parse.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_parse.dir/build: test_parse

.PHONY : CMakeFiles/test_parse.dir/build

CMakeFiles/test_parse.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_parse.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_parse.dir/clean

CMakeFiles/test_parse.dir/depend:
	cd /home/fastF/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fastF /home/fastF /home/fastF/build /home/fastF/build /home/fastF/build/CMakeFiles/test_parse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_parse.dir/depend
