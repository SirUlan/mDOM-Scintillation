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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/m_unla02/Documents/geant4/MyGeant4/6

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/m_unla02/Documents/geant4/MyGeant4/6

# Include any dependencies generated for this target.
include CMakeFiles/mdom_K40.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mdom_K40.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mdom_K40.dir/flags.make

CMakeFiles/mdom_K40.dir/mdom_K40.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/mdom_K40.cc.o: mdom_K40.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/mdom_K40.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/mdom_K40.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/mdom_K40.cc

CMakeFiles/mdom_K40.dir/mdom_K40.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/mdom_K40.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/mdom_K40.cc > CMakeFiles/mdom_K40.dir/mdom_K40.cc.i

CMakeFiles/mdom_K40.dir/mdom_K40.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/mdom_K40.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/mdom_K40.cc -o CMakeFiles/mdom_K40.dir/mdom_K40.cc.s

CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.requires

CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.provides: CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.provides

CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.provides.build: CMakeFiles/mdom_K40.dir/mdom_K40.cc.o

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o: src/mdomRunAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomRunAction.cc

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomRunAction.cc > CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomRunAction.cc -o CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o: src/mdomTrackingAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomTrackingAction.cc

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomTrackingAction.cc > CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomTrackingAction.cc -o CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o: src/SteppingMessenger.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/SteppingMessenger.cc

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/SteppingMessenger.cc > CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.i

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/SteppingMessenger.cc -o CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.s

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.requires

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.provides: CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.provides

CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o: src/mdomDetectorConstruction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomDetectorConstruction.cc

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomDetectorConstruction.cc > CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomDetectorConstruction.cc -o CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o: src/mdomSteppingVerbose.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingVerbose.cc

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingVerbose.cc > CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.i

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingVerbose.cc -o CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.s

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o: src/mdomSteppingAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingAction.cc

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingAction.cc > CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomSteppingAction.cc -o CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o: src/mdomEventAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomEventAction.cc

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomEventAction.cc > CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomEventAction.cc -o CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o: src/mdomAnalysisManager.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomAnalysisManager.cc

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomAnalysisManager.cc > CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.i

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomAnalysisManager.cc -o CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.s

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o: src/mdomPhysicsList.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPhysicsList.cc

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPhysicsList.cc > CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.i

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPhysicsList.cc -o CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.s

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o: CMakeFiles/mdom_K40.dir/flags.make
CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o: src/mdomPrimaryGeneratorAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o -c /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPrimaryGeneratorAction.cc

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPrimaryGeneratorAction.cc > CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.i

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/m_unla02/Documents/geant4/MyGeant4/6/src/mdomPrimaryGeneratorAction.cc -o CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.s

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.requires:
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.requires

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.provides: CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/mdom_K40.dir/build.make CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.provides

CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o

# Object files for target mdom_K40
mdom_K40_OBJECTS = \
"CMakeFiles/mdom_K40.dir/mdom_K40.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o" \
"CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o" \
"CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o"

# External object files for target mdom_K40
mdom_K40_EXTERNAL_OBJECTS =

mdom_K40: CMakeFiles/mdom_K40.dir/mdom_K40.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o
mdom_K40: CMakeFiles/mdom_K40.dir/build.make
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4Tree.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4FR.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4GMocren.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4visHepRep.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4RayTracer.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4VRML.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4OpenGL.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4gl2ps.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4vis_management.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4modeling.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4interfaces.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4persistency.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4analysis.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4error_propagation.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4readout.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4physicslists.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4run.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4event.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4tracking.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4parmodels.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4processes.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4digits_hits.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4track.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4particles.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4geometry.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4materials.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4graphics_reps.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4intercoms.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4global.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4clhep.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4zlib.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4FR.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4vis_management.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4modeling.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libXm.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libSM.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libICE.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libX11.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libXext.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libXmu.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libGLU.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libGL.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libQtGui.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libQtCore.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4run.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4event.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4tracking.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4processes.so
mdom_K40: /usr/lib/x86_64-linux-gnu/libexpat.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4digits_hits.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4track.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4particles.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4geometry.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4materials.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4graphics_reps.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4intercoms.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4global.so
mdom_K40: /iceserver/geant4.9.6.p04_install/lib/libG4clhep.so
mdom_K40: CMakeFiles/mdom_K40.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable mdom_K40"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mdom_K40.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mdom_K40.dir/build: mdom_K40
.PHONY : CMakeFiles/mdom_K40.dir/build

CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/mdom_K40.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomRunAction.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomTrackingAction.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/SteppingMessenger.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomDetectorConstruction.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomSteppingVerbose.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomSteppingAction.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomEventAction.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomAnalysisManager.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomPhysicsList.cc.o.requires
CMakeFiles/mdom_K40.dir/requires: CMakeFiles/mdom_K40.dir/src/mdomPrimaryGeneratorAction.cc.o.requires
.PHONY : CMakeFiles/mdom_K40.dir/requires

CMakeFiles/mdom_K40.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mdom_K40.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mdom_K40.dir/clean

CMakeFiles/mdom_K40.dir/depend:
	cd /home/m_unla02/Documents/geant4/MyGeant4/6 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/m_unla02/Documents/geant4/MyGeant4/6 /home/m_unla02/Documents/geant4/MyGeant4/6 /home/m_unla02/Documents/geant4/MyGeant4/6 /home/m_unla02/Documents/geant4/MyGeant4/6 /home/m_unla02/Documents/geant4/MyGeant4/6/CMakeFiles/mdom_K40.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mdom_K40.dir/depend

