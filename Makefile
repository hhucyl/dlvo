# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/user/workspace1/dlvo/code/dlvo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/workspace1/dlvo/code/dlvo

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/user/workspace1/dlvo/code/dlvo/CMakeFiles /home/user/workspace1/dlvo/code/dlvo/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/user/workspace1/dlvo/code/dlvo/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named test_rw1

# Build rule for target.
test_rw1: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_rw1
.PHONY : test_rw1

# fast build rule for target.
test_rw1/fast:
	$(MAKE) -f CMakeFiles/test_rw1.dir/build.make CMakeFiles/test_rw1.dir/build
.PHONY : test_rw1/fast

#=============================================================================
# Target rules for targets named test_reflect

# Build rule for target.
test_reflect: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_reflect
.PHONY : test_reflect

# fast build rule for target.
test_reflect/fast:
	$(MAKE) -f CMakeFiles/test_reflect.dir/build.make CMakeFiles/test_reflect.dir/build
.PHONY : test_reflect/fast

#=============================================================================
# Target rules for targets named test_swi

# Build rule for target.
test_swi: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_swi
.PHONY : test_swi

# fast build rule for target.
test_swi/fast:
	$(MAKE) -f CMakeFiles/test_swi.dir/build.make CMakeFiles/test_swi.dir/build
.PHONY : test_swi/fast

#=============================================================================
# Target rules for targets named test_2

# Build rule for target.
test_2: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_2
.PHONY : test_2

# fast build rule for target.
test_2/fast:
	$(MAKE) -f CMakeFiles/test_2.dir/build.make CMakeFiles/test_2.dir/build
.PHONY : test_2/fast

#=============================================================================
# Target rules for targets named test_cong

# Build rule for target.
test_cong: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_cong
.PHONY : test_cong

# fast build rule for target.
test_cong/fast:
	$(MAKE) -f CMakeFiles/test_cong.dir/build.make CMakeFiles/test_cong.dir/build
.PHONY : test_cong/fast

#=============================================================================
# Target rules for targets named test_fine

# Build rule for target.
test_fine: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_fine
.PHONY : test_fine

# fast build rule for target.
test_fine/fast:
	$(MAKE) -f CMakeFiles/test_fine.dir/build.make CMakeFiles/test_fine.dir/build
.PHONY : test_fine/fast

#=============================================================================
# Target rules for targets named test_swi_h

# Build rule for target.
test_swi_h: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_swi_h
.PHONY : test_swi_h

# fast build rule for target.
test_swi_h/fast:
	$(MAKE) -f CMakeFiles/test_swi_h.dir/build.make CMakeFiles/test_swi_h.dir/build
.PHONY : test_swi_h/fast

#=============================================================================
# Target rules for targets named test_wu

# Build rule for target.
test_wu: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_wu
.PHONY : test_wu

# fast build rule for target.
test_wu/fast:
	$(MAKE) -f CMakeFiles/test_wu.dir/build.make CMakeFiles/test_wu.dir/build
.PHONY : test_wu/fast

#=============================================================================
# Target rules for targets named test_wu2

# Build rule for target.
test_wu2: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_wu2
.PHONY : test_wu2

# fast build rule for target.
test_wu2/fast:
	$(MAKE) -f CMakeFiles/test_wu2.dir/build.make CMakeFiles/test_wu2.dir/build
.PHONY : test_wu2/fast

#=============================================================================
# Target rules for targets named test_swi_rw

# Build rule for target.
test_swi_rw: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_swi_rw
.PHONY : test_swi_rw

# fast build rule for target.
test_swi_rw/fast:
	$(MAKE) -f CMakeFiles/test_swi_rw.dir/build.make CMakeFiles/test_swi_rw.dir/build
.PHONY : test_swi_rw/fast

#=============================================================================
# Target rules for targets named test_rw

# Build rule for target.
test_rw: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_rw
.PHONY : test_rw

# fast build rule for target.
test_rw/fast:
	$(MAKE) -f CMakeFiles/test_rw.dir/build.make CMakeFiles/test_rw.dir/build
.PHONY : test_rw/fast

#=============================================================================
# Target rules for targets named test_wu3

# Build rule for target.
test_wu3: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_wu3
.PHONY : test_wu3

# fast build rule for target.
test_wu3/fast:
	$(MAKE) -f CMakeFiles/test_wu3.dir/build.make CMakeFiles/test_wu3.dir/build
.PHONY : test_wu3/fast

test_2.o: test_2.cpp.o

.PHONY : test_2.o

# target to build an object file
test_2.cpp.o:
	$(MAKE) -f CMakeFiles/test_2.dir/build.make CMakeFiles/test_2.dir/test_2.cpp.o
.PHONY : test_2.cpp.o

test_2.i: test_2.cpp.i

.PHONY : test_2.i

# target to preprocess a source file
test_2.cpp.i:
	$(MAKE) -f CMakeFiles/test_2.dir/build.make CMakeFiles/test_2.dir/test_2.cpp.i
.PHONY : test_2.cpp.i

test_2.s: test_2.cpp.s

.PHONY : test_2.s

# target to generate assembly for a file
test_2.cpp.s:
	$(MAKE) -f CMakeFiles/test_2.dir/build.make CMakeFiles/test_2.dir/test_2.cpp.s
.PHONY : test_2.cpp.s

test_cong.o: test_cong.cpp.o

.PHONY : test_cong.o

# target to build an object file
test_cong.cpp.o:
	$(MAKE) -f CMakeFiles/test_cong.dir/build.make CMakeFiles/test_cong.dir/test_cong.cpp.o
.PHONY : test_cong.cpp.o

test_cong.i: test_cong.cpp.i

.PHONY : test_cong.i

# target to preprocess a source file
test_cong.cpp.i:
	$(MAKE) -f CMakeFiles/test_cong.dir/build.make CMakeFiles/test_cong.dir/test_cong.cpp.i
.PHONY : test_cong.cpp.i

test_cong.s: test_cong.cpp.s

.PHONY : test_cong.s

# target to generate assembly for a file
test_cong.cpp.s:
	$(MAKE) -f CMakeFiles/test_cong.dir/build.make CMakeFiles/test_cong.dir/test_cong.cpp.s
.PHONY : test_cong.cpp.s

test_fine.o: test_fine.cpp.o

.PHONY : test_fine.o

# target to build an object file
test_fine.cpp.o:
	$(MAKE) -f CMakeFiles/test_fine.dir/build.make CMakeFiles/test_fine.dir/test_fine.cpp.o
.PHONY : test_fine.cpp.o

test_fine.i: test_fine.cpp.i

.PHONY : test_fine.i

# target to preprocess a source file
test_fine.cpp.i:
	$(MAKE) -f CMakeFiles/test_fine.dir/build.make CMakeFiles/test_fine.dir/test_fine.cpp.i
.PHONY : test_fine.cpp.i

test_fine.s: test_fine.cpp.s

.PHONY : test_fine.s

# target to generate assembly for a file
test_fine.cpp.s:
	$(MAKE) -f CMakeFiles/test_fine.dir/build.make CMakeFiles/test_fine.dir/test_fine.cpp.s
.PHONY : test_fine.cpp.s

test_reflect.o: test_reflect.cpp.o

.PHONY : test_reflect.o

# target to build an object file
test_reflect.cpp.o:
	$(MAKE) -f CMakeFiles/test_reflect.dir/build.make CMakeFiles/test_reflect.dir/test_reflect.cpp.o
.PHONY : test_reflect.cpp.o

test_reflect.i: test_reflect.cpp.i

.PHONY : test_reflect.i

# target to preprocess a source file
test_reflect.cpp.i:
	$(MAKE) -f CMakeFiles/test_reflect.dir/build.make CMakeFiles/test_reflect.dir/test_reflect.cpp.i
.PHONY : test_reflect.cpp.i

test_reflect.s: test_reflect.cpp.s

.PHONY : test_reflect.s

# target to generate assembly for a file
test_reflect.cpp.s:
	$(MAKE) -f CMakeFiles/test_reflect.dir/build.make CMakeFiles/test_reflect.dir/test_reflect.cpp.s
.PHONY : test_reflect.cpp.s

test_rw.o: test_rw.cpp.o

.PHONY : test_rw.o

# target to build an object file
test_rw.cpp.o:
	$(MAKE) -f CMakeFiles/test_rw.dir/build.make CMakeFiles/test_rw.dir/test_rw.cpp.o
.PHONY : test_rw.cpp.o

test_rw.i: test_rw.cpp.i

.PHONY : test_rw.i

# target to preprocess a source file
test_rw.cpp.i:
	$(MAKE) -f CMakeFiles/test_rw.dir/build.make CMakeFiles/test_rw.dir/test_rw.cpp.i
.PHONY : test_rw.cpp.i

test_rw.s: test_rw.cpp.s

.PHONY : test_rw.s

# target to generate assembly for a file
test_rw.cpp.s:
	$(MAKE) -f CMakeFiles/test_rw.dir/build.make CMakeFiles/test_rw.dir/test_rw.cpp.s
.PHONY : test_rw.cpp.s

test_rw1.o: test_rw1.cpp.o

.PHONY : test_rw1.o

# target to build an object file
test_rw1.cpp.o:
	$(MAKE) -f CMakeFiles/test_rw1.dir/build.make CMakeFiles/test_rw1.dir/test_rw1.cpp.o
.PHONY : test_rw1.cpp.o

test_rw1.i: test_rw1.cpp.i

.PHONY : test_rw1.i

# target to preprocess a source file
test_rw1.cpp.i:
	$(MAKE) -f CMakeFiles/test_rw1.dir/build.make CMakeFiles/test_rw1.dir/test_rw1.cpp.i
.PHONY : test_rw1.cpp.i

test_rw1.s: test_rw1.cpp.s

.PHONY : test_rw1.s

# target to generate assembly for a file
test_rw1.cpp.s:
	$(MAKE) -f CMakeFiles/test_rw1.dir/build.make CMakeFiles/test_rw1.dir/test_rw1.cpp.s
.PHONY : test_rw1.cpp.s

test_swi.o: test_swi.cpp.o

.PHONY : test_swi.o

# target to build an object file
test_swi.cpp.o:
	$(MAKE) -f CMakeFiles/test_swi.dir/build.make CMakeFiles/test_swi.dir/test_swi.cpp.o
.PHONY : test_swi.cpp.o

test_swi.i: test_swi.cpp.i

.PHONY : test_swi.i

# target to preprocess a source file
test_swi.cpp.i:
	$(MAKE) -f CMakeFiles/test_swi.dir/build.make CMakeFiles/test_swi.dir/test_swi.cpp.i
.PHONY : test_swi.cpp.i

test_swi.s: test_swi.cpp.s

.PHONY : test_swi.s

# target to generate assembly for a file
test_swi.cpp.s:
	$(MAKE) -f CMakeFiles/test_swi.dir/build.make CMakeFiles/test_swi.dir/test_swi.cpp.s
.PHONY : test_swi.cpp.s

test_swi_h.o: test_swi_h.cpp.o

.PHONY : test_swi_h.o

# target to build an object file
test_swi_h.cpp.o:
	$(MAKE) -f CMakeFiles/test_swi_h.dir/build.make CMakeFiles/test_swi_h.dir/test_swi_h.cpp.o
.PHONY : test_swi_h.cpp.o

test_swi_h.i: test_swi_h.cpp.i

.PHONY : test_swi_h.i

# target to preprocess a source file
test_swi_h.cpp.i:
	$(MAKE) -f CMakeFiles/test_swi_h.dir/build.make CMakeFiles/test_swi_h.dir/test_swi_h.cpp.i
.PHONY : test_swi_h.cpp.i

test_swi_h.s: test_swi_h.cpp.s

.PHONY : test_swi_h.s

# target to generate assembly for a file
test_swi_h.cpp.s:
	$(MAKE) -f CMakeFiles/test_swi_h.dir/build.make CMakeFiles/test_swi_h.dir/test_swi_h.cpp.s
.PHONY : test_swi_h.cpp.s

test_swi_rw.o: test_swi_rw.cpp.o

.PHONY : test_swi_rw.o

# target to build an object file
test_swi_rw.cpp.o:
	$(MAKE) -f CMakeFiles/test_swi_rw.dir/build.make CMakeFiles/test_swi_rw.dir/test_swi_rw.cpp.o
.PHONY : test_swi_rw.cpp.o

test_swi_rw.i: test_swi_rw.cpp.i

.PHONY : test_swi_rw.i

# target to preprocess a source file
test_swi_rw.cpp.i:
	$(MAKE) -f CMakeFiles/test_swi_rw.dir/build.make CMakeFiles/test_swi_rw.dir/test_swi_rw.cpp.i
.PHONY : test_swi_rw.cpp.i

test_swi_rw.s: test_swi_rw.cpp.s

.PHONY : test_swi_rw.s

# target to generate assembly for a file
test_swi_rw.cpp.s:
	$(MAKE) -f CMakeFiles/test_swi_rw.dir/build.make CMakeFiles/test_swi_rw.dir/test_swi_rw.cpp.s
.PHONY : test_swi_rw.cpp.s

test_wu.o: test_wu.cpp.o

.PHONY : test_wu.o

# target to build an object file
test_wu.cpp.o:
	$(MAKE) -f CMakeFiles/test_wu.dir/build.make CMakeFiles/test_wu.dir/test_wu.cpp.o
.PHONY : test_wu.cpp.o

test_wu.i: test_wu.cpp.i

.PHONY : test_wu.i

# target to preprocess a source file
test_wu.cpp.i:
	$(MAKE) -f CMakeFiles/test_wu.dir/build.make CMakeFiles/test_wu.dir/test_wu.cpp.i
.PHONY : test_wu.cpp.i

test_wu.s: test_wu.cpp.s

.PHONY : test_wu.s

# target to generate assembly for a file
test_wu.cpp.s:
	$(MAKE) -f CMakeFiles/test_wu.dir/build.make CMakeFiles/test_wu.dir/test_wu.cpp.s
.PHONY : test_wu.cpp.s

test_wu2.o: test_wu2.cpp.o

.PHONY : test_wu2.o

# target to build an object file
test_wu2.cpp.o:
	$(MAKE) -f CMakeFiles/test_wu2.dir/build.make CMakeFiles/test_wu2.dir/test_wu2.cpp.o
.PHONY : test_wu2.cpp.o

test_wu2.i: test_wu2.cpp.i

.PHONY : test_wu2.i

# target to preprocess a source file
test_wu2.cpp.i:
	$(MAKE) -f CMakeFiles/test_wu2.dir/build.make CMakeFiles/test_wu2.dir/test_wu2.cpp.i
.PHONY : test_wu2.cpp.i

test_wu2.s: test_wu2.cpp.s

.PHONY : test_wu2.s

# target to generate assembly for a file
test_wu2.cpp.s:
	$(MAKE) -f CMakeFiles/test_wu2.dir/build.make CMakeFiles/test_wu2.dir/test_wu2.cpp.s
.PHONY : test_wu2.cpp.s

test_wu3.o: test_wu3.cpp.o

.PHONY : test_wu3.o

# target to build an object file
test_wu3.cpp.o:
	$(MAKE) -f CMakeFiles/test_wu3.dir/build.make CMakeFiles/test_wu3.dir/test_wu3.cpp.o
.PHONY : test_wu3.cpp.o

test_wu3.i: test_wu3.cpp.i

.PHONY : test_wu3.i

# target to preprocess a source file
test_wu3.cpp.i:
	$(MAKE) -f CMakeFiles/test_wu3.dir/build.make CMakeFiles/test_wu3.dir/test_wu3.cpp.i
.PHONY : test_wu3.cpp.i

test_wu3.s: test_wu3.cpp.s

.PHONY : test_wu3.s

# target to generate assembly for a file
test_wu3.cpp.s:
	$(MAKE) -f CMakeFiles/test_wu3.dir/build.make CMakeFiles/test_wu3.dir/test_wu3.cpp.s
.PHONY : test_wu3.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... test_rw1"
	@echo "... test_reflect"
	@echo "... test_swi"
	@echo "... test_2"
	@echo "... test_cong"
	@echo "... test_fine"
	@echo "... test_swi_h"
	@echo "... test_wu"
	@echo "... test_wu2"
	@echo "... test_swi_rw"
	@echo "... test_rw"
	@echo "... test_wu3"
	@echo "... test_2.o"
	@echo "... test_2.i"
	@echo "... test_2.s"
	@echo "... test_cong.o"
	@echo "... test_cong.i"
	@echo "... test_cong.s"
	@echo "... test_fine.o"
	@echo "... test_fine.i"
	@echo "... test_fine.s"
	@echo "... test_reflect.o"
	@echo "... test_reflect.i"
	@echo "... test_reflect.s"
	@echo "... test_rw.o"
	@echo "... test_rw.i"
	@echo "... test_rw.s"
	@echo "... test_rw1.o"
	@echo "... test_rw1.i"
	@echo "... test_rw1.s"
	@echo "... test_swi.o"
	@echo "... test_swi.i"
	@echo "... test_swi.s"
	@echo "... test_swi_h.o"
	@echo "... test_swi_h.i"
	@echo "... test_swi_h.s"
	@echo "... test_swi_rw.o"
	@echo "... test_swi_rw.i"
	@echo "... test_swi_rw.s"
	@echo "... test_wu.o"
	@echo "... test_wu.i"
	@echo "... test_wu.s"
	@echo "... test_wu2.o"
	@echo "... test_wu2.i"
	@echo "... test_wu2.s"
	@echo "... test_wu3.o"
	@echo "... test_wu3.i"
	@echo "... test_wu3.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

