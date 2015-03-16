# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake

# The command to remove a file.
RM = /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratch2/chongchu/bamtools-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch2/chongchu/bamtools-master/build

# Utility rule file for APIHeaders.

# Include the progress variables for this target.
include src/api/CMakeFiles/APIHeaders.dir/progress.make

src/api/CMakeFiles/APIHeaders:
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch2/chongchu/bamtools-master/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Exporting APIHeaders"

APIHeaders: src/api/CMakeFiles/APIHeaders
APIHeaders: src/api/CMakeFiles/APIHeaders.dir/build.make
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/api_global.h /scratch2/chongchu/bamtools-master/include/api/api_global.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamAlgorithms.h /scratch2/chongchu/bamtools-master/include/api/BamAlgorithms.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamAlignment.h /scratch2/chongchu/bamtools-master/include/api/BamAlignment.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamAux.h /scratch2/chongchu/bamtools-master/include/api/BamAux.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamConstants.h /scratch2/chongchu/bamtools-master/include/api/BamConstants.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamIndex.h /scratch2/chongchu/bamtools-master/include/api/BamIndex.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamMultiReader.h /scratch2/chongchu/bamtools-master/include/api/BamMultiReader.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamReader.h /scratch2/chongchu/bamtools-master/include/api/BamReader.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/BamWriter.h /scratch2/chongchu/bamtools-master/include/api/BamWriter.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/IBamIODevice.h /scratch2/chongchu/bamtools-master/include/api/IBamIODevice.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamConstants.h /scratch2/chongchu/bamtools-master/include/api/SamConstants.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamHeader.h /scratch2/chongchu/bamtools-master/include/api/SamHeader.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamProgram.h /scratch2/chongchu/bamtools-master/include/api/SamProgram.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamProgramChain.h /scratch2/chongchu/bamtools-master/include/api/SamProgramChain.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamReadGroup.h /scratch2/chongchu/bamtools-master/include/api/SamReadGroup.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamReadGroupDictionary.h /scratch2/chongchu/bamtools-master/include/api/SamReadGroupDictionary.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamSequence.h /scratch2/chongchu/bamtools-master/include/api/SamSequence.h
	cd /scratch2/chongchu/bamtools-master/build/src/api && /scratch2/chongchu/cmake-2.8.10.2-Linux-i386/bin/cmake -E copy_if_different /scratch2/chongchu/bamtools-master/src/api/SamSequenceDictionary.h /scratch2/chongchu/bamtools-master/include/api/SamSequenceDictionary.h
.PHONY : APIHeaders

# Rule to build all files generated by this target.
src/api/CMakeFiles/APIHeaders.dir/build: APIHeaders
.PHONY : src/api/CMakeFiles/APIHeaders.dir/build

src/api/CMakeFiles/APIHeaders.dir/clean:
	cd /scratch2/chongchu/bamtools-master/build/src/api && $(CMAKE_COMMAND) -P CMakeFiles/APIHeaders.dir/cmake_clean.cmake
.PHONY : src/api/CMakeFiles/APIHeaders.dir/clean

src/api/CMakeFiles/APIHeaders.dir/depend:
	cd /scratch2/chongchu/bamtools-master/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch2/chongchu/bamtools-master /scratch2/chongchu/bamtools-master/src/api /scratch2/chongchu/bamtools-master/build /scratch2/chongchu/bamtools-master/build/src/api /scratch2/chongchu/bamtools-master/build/src/api/CMakeFiles/APIHeaders.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/api/CMakeFiles/APIHeaders.dir/depend

