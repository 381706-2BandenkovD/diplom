cmake_minimum_required(VERSION 3.7)

project( Sience_Lab )
FIND_PACKAGE(PythonInterp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(3rdparty/gtest)
add_definitions(-DDATA_FOLDER="${PROJECT_SOURCE_DIR}/result")

if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python)
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  else()	
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  endif()
else()
    message("Python not found")
endif()

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file(GLOB modules "modules/*")
foreach(module ${modules})
    get_filename_component(module_name ${module} NAME)
    string(REGEX REPLACE "[0-9]+_" "" module_name ${module_name})

    # Build library
    file(GLOB SOURCES "${module}/src/*.cpp")
    if (SOURCES)
        add_library(${module_name} STATIC ${SOURCES})
        target_include_directories(${module_name} PUBLIC "${module}/include" ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
        target_link_libraries(${module_name} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    endif()
    
    # Compile test
    file(GLOB TEST_SOURCES "${module}/*.cpp")
    if (TEST_SOURCES)
        PYTHON_ADD_MODULE(${module_name} "${TEST_SOURCES}")
        target_include_directories(${module_name} PUBLIC ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
        target_link_libraries(${module_name} PRIVATE math_model PUBLIC ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    endif()

    file(GLOB PYTEST_FILE "${module}/test_${module_name}.py")
    if (PYTEST_FILE)
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E copy
                ${PYTEST_FILE}
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${module_name})
    endif()
endforeach()


