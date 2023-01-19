cmake_minimum_required(VERSION 3.18)
project(cpp-d4)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)
  	find_package(LAPACKE REQUIRED)
	find_package(CBLAS CONFIG REQUIRED)
endif(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)


IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -pedantic \
                                           -Wextra -Wcast-align -Wcast-qual  -Wchar-subscripts  \
                                           -Wcomment -Wdisabled-optimization \
                                           -Wformat  -Wformat=2 -Wformat-nonliteral -Wformat-security\
                                           -Wformat-y2k  -Wimport  -Winit-self  -Winline -Winvalid-pch\
                                           -Wunsafe-loop-optimizations  -Wmissing-braces\
                                           -Wmissing-field-initializers -Wmissing-format-attribute  \
                                           -Wmissing-include-dirs -Wmissing-noreturn -Wpacked  -Wparentheses\
                                           -Wpointer-arith -Wredundant-decls -Wsequence-point\
                                           -Wsign-compare  -Wstack-protector -Wstrict-aliasing\
                                           -Wstrict-aliasing=2 -Wswitch -Wsuggest-override\
                                           -Wtrigraphs  -Wuninitialized -Wunknown-pragmas  -Wunreachable-code\
                                           -Wunused -Wunused-function  -Wunused-label  -Wunused-parameter\
                                           -Wunused-value  -Wunused-variable  -Wvariadic-macros\
                                           -Wvolatile-register-var  -Wwrite-strings -Wdeprecated-declarations\
                                           -Wno-error=unused-local-typedefs -Wno-error=enum-compare -Wno-narrowing -Werror=return-type -g")
set (GCC ON)
if(WIN32) # Check if we are on Windows
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
endif(WIN32)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR} )
if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
	include_directories(${D4INCLUDE_DIR})
endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)


set(cpp_d4_SRC
  src/damping/atm.cpp
  src/damping/rational.cpp
  src/dftd_cutoff.cpp
  src/dftd_damping.cpp
  src/dftd_dispersion.cpp
  src/dftd_eeq.cpp 
  src/dftd_model.cpp
  src/dftd_ncoord.cpp
  src/dftd_readxyz.cpp
  )
  
add_library(libcpp_d4  ${cpp_d4_SRC})

add_executable(cpp-d4
  src/program_dftd.cpp
)

add_executable(cpp-d4_test
  test/main.cpp
  test/test_grad.cpp
  test/molecules.cpp
  test/test_disp.cpp
  test/test_ncoord.cpp
  test/test_param.cpp
  test/util.cpp
)

if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
	link_directories(${D4LIBS_DIR})

	target_link_libraries(cpp-d4 libcpp_d4 cblas  lapacke pthread)
	target_link_libraries(cpp-d4_test libcpp_d4 cblas lapacke pthread)
else()
	target_link_libraries(cpp-d4 libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} pthread)
	target_link_libraries(cpp-d4_test libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} pthread)
endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)

include(CTest)
add_test(NAME disp COMMAND cpp-d4_test disp)
add_test(NAME grad COMMAND cpp-d4_test grad)
add_test(NAME ncoord COMMAND cpp-d4_test ncoord)
add_test(NAME param COMMAND cpp-d4_test param)

if(WIN32) # Check if we are on Windows
else()
     target_link_libraries(cpp-d4 dl )
endif(WIN32)
