
# Se compila todo el proyecto usando gfortran
cmake_minimum_required(VERSION 3.16)

project(mcfort)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -Ofast -march=native -funroll-loops")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compilar
#
# Esta linea junta a todos los archivos y al programa principal
file(GLOB_RECURSE sources include/*.f90 src/*.f90)
add_executable(mcfort ${sources})