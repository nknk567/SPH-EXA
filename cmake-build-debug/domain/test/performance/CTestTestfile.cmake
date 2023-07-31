# CMake generated Testfile for 
# Source directory: /Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/test/performance
# Build directory: /Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/cmake-build-debug/domain/test/performance
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(scan_perf "/opt/homebrew/bin/mpiexec" "-n" "1" "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/cmake-build-debug/domain/test/performance/scan_perf")
set_tests_properties(scan_perf PROPERTIES  _BACKTRACE_TRIPLES "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;8;cstone_add_test;/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;0;")
