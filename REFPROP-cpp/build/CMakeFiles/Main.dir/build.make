# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = X:\CoBraPython\REFPROP-cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = X:\CoBraPython\REFPROP-cpp\build

# Include any dependencies generated for this target.
include CMakeFiles/Main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Main.dir/flags.make

CMakeFiles/Main.dir/main.cpp.obj: CMakeFiles/Main.dir/flags.make
CMakeFiles/Main.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=X:\CoBraPython\REFPROP-cpp\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Main.dir/main.cpp.obj"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Main.dir\main.cpp.obj -c X:\CoBraPython\REFPROP-cpp\main.cpp

CMakeFiles/Main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Main.dir/main.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E X:\CoBraPython\REFPROP-cpp\main.cpp > CMakeFiles\Main.dir\main.cpp.i

CMakeFiles/Main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Main.dir/main.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S X:\CoBraPython\REFPROP-cpp\main.cpp -o CMakeFiles\Main.dir\main.cpp.s

CMakeFiles/Main.dir/other.cpp.obj: CMakeFiles/Main.dir/flags.make
CMakeFiles/Main.dir/other.cpp.obj: ../other.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=X:\CoBraPython\REFPROP-cpp\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Main.dir/other.cpp.obj"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Main.dir\other.cpp.obj -c X:\CoBraPython\REFPROP-cpp\other.cpp

CMakeFiles/Main.dir/other.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Main.dir/other.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E X:\CoBraPython\REFPROP-cpp\other.cpp > CMakeFiles\Main.dir\other.cpp.i

CMakeFiles/Main.dir/other.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Main.dir/other.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S X:\CoBraPython\REFPROP-cpp\other.cpp -o CMakeFiles\Main.dir\other.cpp.s

# Object files for target Main
Main_OBJECTS = \
"CMakeFiles/Main.dir/main.cpp.obj" \
"CMakeFiles/Main.dir/other.cpp.obj"

# External object files for target Main
Main_EXTERNAL_OBJECTS =

Main.exe: CMakeFiles/Main.dir/main.cpp.obj
Main.exe: CMakeFiles/Main.dir/other.cpp.obj
Main.exe: CMakeFiles/Main.dir/build.make
Main.exe: CMakeFiles/Main.dir/linklibs.rsp
Main.exe: CMakeFiles/Main.dir/objects1.rsp
Main.exe: CMakeFiles/Main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=X:\CoBraPython\REFPROP-cpp\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Main.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Main.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Main.dir/build: Main.exe

.PHONY : CMakeFiles/Main.dir/build

CMakeFiles/Main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Main.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Main.dir/clean

CMakeFiles/Main.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" X:\CoBraPython\REFPROP-cpp X:\CoBraPython\REFPROP-cpp X:\CoBraPython\REFPROP-cpp\build X:\CoBraPython\REFPROP-cpp\build X:\CoBraPython\REFPROP-cpp\build\CMakeFiles\Main.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Main.dir/depend

