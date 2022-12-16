cmake_minimum_required(VERSION 3.4)

project(fesom Fortran)

option(ENABLE_OPENACC "compile with OpenACC support" ON)
set(GPU_COMPUTE_CAPABILITY "cc80" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")
set(GPU_FLAGS "cuda11.7,${GPU_COMPUTE_CAPABILITY}" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")

option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_Fortran ${src_home}/*.F90)

#add_library(${PROJECT_NAME}_C ${sources_C})
#target_compile_definitions(${PROJECT_NAME}_C PRIVATE PARMS USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
#target_link_libraries(${PROJECT_NAME}_C parms) #metis

# create our binary (set its name to name of this project)
add_executable(${PROJECT_NAME} ${sources_Fortran})
if(${VERBOSE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
        target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
        target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
                target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
        endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC)
   if(${ENABLE_OPENACC})
       set(GPU_COMPUTE_CAPABILITY "cc80" CACHE STRING "GPU Arch for PGI Compiler (cc35,cc50,cc60,cc70,cc80,...)")
       set(GPU_FLAGS "cuda11.7,${GPU_COMPUTE_CAPABILITY}" CACHE STRING "Flags to be passed to \"-gpu\" when building gpu code")
       target_compile_options(${PROJECT_NAME} PRIVATE -Mallocatable=95 -Mr8 -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=${GPU_FLAGS})
       # target_compile_options(${PROJECT_NAME} PRIVATE -Mallocatable=95 -Mr8 -pgf90libs -Minfo=all -O3 -fast -acc=verystrict -gpu=${GPU_FLAGS})
       set(CMAKE_EXE_LINKER_FLAGS "-acc=verystrict -gpu=${GPU_FLAGS}")
   else()
       target_compile_options(${PROJECT_NAME} PRIVATE -Mallocatable=95 -Mr8 -pgf90libs -Mnofma)
   endif()
endif()

target_link_libraries(${PROJECT_NAME} nvhpcwrapnvtx)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fdiff" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
