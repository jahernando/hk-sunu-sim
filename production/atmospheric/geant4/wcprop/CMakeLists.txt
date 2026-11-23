cmake_minimum_required(VERSION 3.16)
project(wcprop CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
add_executable(wcprop wcprop.cc)
target_link_libraries(wcprop ${Geant4_LIBRARIES})
if(APPLE)
  # In conda environments the lib dir can end up twice in LC_RPATH, which
  # recent macOS dyld rejects. Remove one copy (no-op if there is no duplicate).
  add_custom_command(TARGET wcprop POST_BUILD
    COMMAND sh -c "n=$(otool -l $<TARGET_FILE:wcprop> | grep -c 'path $ENV{CONDA_PREFIX}/lib '); if [ \"$n\" -gt 1 ]; then install_name_tool -delete_rpath $ENV{CONDA_PREFIX}/lib $<TARGET_FILE:wcprop>; fi"
    VERBATIM)
endif()
