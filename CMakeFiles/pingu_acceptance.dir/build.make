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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ray/Programme/geant4/pingu.exp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ray/Programme/geant4/pingu.exp

# Include any dependencies generated for this target.
include CMakeFiles/pingu_acceptance.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pingu_acceptance.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pingu_acceptance.dir/flags.make

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o: pingu_acceptance.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o -c /home/ray/Programme/geant4/pingu.exp/pingu_acceptance.cc

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/pingu_acceptance.cc > CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.i

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/pingu_acceptance.cc -o CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.s

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.requires

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.provides: CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.provides

CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o: src/pinguTrackingAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguTrackingAction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguTrackingAction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguTrackingAction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o: src/pinguRunAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguRunAction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguRunAction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguRunAction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o: src/pinguSteppingVerbose.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingVerbose.cc

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingVerbose.cc > CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingVerbose.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o: src/pinguEventAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguEventAction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguEventAction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguEventAction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o: src/pinguPrimaryGeneratorAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguPrimaryGeneratorAction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguPrimaryGeneratorAction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguPrimaryGeneratorAction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o: src/pinguSteppingAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingAction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingAction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguSteppingAction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o: src/pinguDetectorConstruction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguDetectorConstruction.cc

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguDetectorConstruction.cc > CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguDetectorConstruction.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o: CMakeFiles/pingu_acceptance.dir/flags.make
CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o: src/pinguPhysicsList.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ray/Programme/geant4/pingu.exp/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o -c /home/ray/Programme/geant4/pingu.exp/src/pinguPhysicsList.cc

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ray/Programme/geant4/pingu.exp/src/pinguPhysicsList.cc > CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.i

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ray/Programme/geant4/pingu.exp/src/pinguPhysicsList.cc -o CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.s

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.requires:
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.requires

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.provides: CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/pingu_acceptance.dir/build.make CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.provides

CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.provides.build: CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o

# Object files for target pingu_acceptance
pingu_acceptance_OBJECTS = \
"CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o" \
"CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o"

# External object files for target pingu_acceptance
pingu_acceptance_EXTERNAL_OBJECTS =

pingu_acceptance: CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4Tree.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4FR.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4GMocren.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4visHepRep.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4RayTracer.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4VRML.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4OpenGL.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4gl2ps.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4vis_management.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4modeling.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4interfaces.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4persistency.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4analysis.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4error_propagation.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4readout.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4physicslists.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4run.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4event.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4tracking.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4parmodels.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4processes.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4digits_hits.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4track.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4particles.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4geometry.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4materials.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4graphics_reps.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4intercoms.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4global.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4clhep.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4zlib.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4FR.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4vis_management.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4modeling.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libGLU.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libGL.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libSM.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libICE.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libX11.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libXext.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libXmu.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4run.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4event.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4tracking.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4processes.so
pingu_acceptance: /usr/lib/i386-linux-gnu/libexpat.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4digits_hits.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4track.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4particles.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4geometry.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4materials.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4graphics_reps.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4intercoms.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4global.so
pingu_acceptance: /home/ray/Programme/geant4/geant4_9_6_p02/lib/libG4clhep.so
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/build.make
pingu_acceptance: CMakeFiles/pingu_acceptance.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable pingu_acceptance"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pingu_acceptance.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pingu_acceptance.dir/build: pingu_acceptance
.PHONY : CMakeFiles/pingu_acceptance.dir/build

CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/pingu_acceptance.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguTrackingAction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguRunAction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingVerbose.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguEventAction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguPrimaryGeneratorAction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguSteppingAction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguDetectorConstruction.cc.o.requires
CMakeFiles/pingu_acceptance.dir/requires: CMakeFiles/pingu_acceptance.dir/src/pinguPhysicsList.cc.o.requires
.PHONY : CMakeFiles/pingu_acceptance.dir/requires

CMakeFiles/pingu_acceptance.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pingu_acceptance.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pingu_acceptance.dir/clean

CMakeFiles/pingu_acceptance.dir/depend:
	cd /home/ray/Programme/geant4/pingu.exp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ray/Programme/geant4/pingu.exp /home/ray/Programme/geant4/pingu.exp /home/ray/Programme/geant4/pingu.exp /home/ray/Programme/geant4/pingu.exp /home/ray/Programme/geant4/pingu.exp/CMakeFiles/pingu_acceptance.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pingu_acceptance.dir/depend

