# CMake generated Testfile for 
# Source directory: /Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/main/test/mpi
# Build directory: /Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/cmake-build-debug/main/test/mpi
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[HDF5IO]=] "/opt/homebrew/bin/mpiexec" "-n" "2" "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/cmake-build-debug/main/test/mpi/hdf5io")
set_tests_properties([=[HDF5IO]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/test/integration_mpi/CMakeLists.txt;10;cstone_add_test;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;7;addMpiTest;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;12;addFrontendMpiTest;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;0;")
