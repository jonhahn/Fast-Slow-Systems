cmake_minimum_required(VERSION 2.8)

# macro for enabling c99 pulled from here: http://stackoverflow.com/a/30564223
macro(use_c99)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")
    endif ()
  else ()
    set (CMAKE_C_STANDARD 99)
  endif ()
endmacro(use_c99)

use_c99()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

set(SOURCE_FILES 

        src/core/rk45.c
        inc/core/rk45.h

        src/core/outputJson.c
        inc/core/outputJson.h

        src/core/parson/parson.c
        inc/core/parson/parson.h

        inc/core/systems.h
        src/core/systems.c

        inc/newtonsMethodforQuadraticSystem.h
        src/variationalPoincare/newtonsMethodforQuadraticSystem.c

        )

include_directories(
        inc
)

add_library(library STATIC ${SOURCE_FILES})

add_executable(testVanDerPol src/tests/testVanDerPol.c)
target_link_libraries(testVanDerPol library m)

add_executable(vanDerPolPoincare src/VDPpoincare/vanDerPolPoincare.c)
target_link_libraries(vanDerPolPoincare library m)

add_executable(testQuadraticPoincare src/variationalPoincare/quadratic.c)
target_link_libraries(testQuadraticPoincare library m)

add_executable(newtonsMethodForQuadraticSystemMain src/variationalPoincare/newtonsMethodForQuadraticSystemMain.c)
target_link_libraries(newtonsMethodForQuadraticSystemMain library m)

add_executable(basinOfAttraction src/variationalPoincare/basinOfAttraction.c)
target_link_libraries(basinOfAttraction library m)

add_executable(onePeriod src/variationalPoincare/onePeriod.c)
target_link_libraries(onePeriod library m)

add_executable(makeBifurcationDiagramVDP src/VDPpoincare/makeBifurcationDiagram.c)
target_link_libraries(makeBifurcationDiagramVDP library m)
