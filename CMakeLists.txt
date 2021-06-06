cmake_minimum_required( VERSION 3.20 )
project( qsimplex VERSION 0.0.1 )

enable_language( CXX )
set( CMAKE_CONFIGURATION_TYPES "Debug;Release;Native" )

set( CMAKE_CXX_FLAGS         "-std=c++20" )
set( CMAKE_CXX_FLAGS_DEBUG   "-Wall -Wextra -Wpedantic -Og -march=core2 -mtune=native -ggdb3" )
set( CMAKE_CXX_FLAGS_RELEASE "-w -O2 -march=core2 -mtune=generic -DNDEBUG" )
set( CMAKE_CXX_FLAGS_NATIVE  "-w -Ofast -march=native -mtune=native -DNDEBUG" )

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release )
endif()
message( "Using build configuration: ${CMAKE_BUILD_TYPE}" )

find_package(Qt6 COMPONENTS Core Gui REQUIRED)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_SOURCE_DIR}/src/ )

file( GLOB_RECURSE MAIN_FILES 
	"${CMAKE_SOURCE_DIR}/src/*.cc"
)

add_executable( qsimplex ${MAIN_FILES} )
target_link_libraries( qsimplex PUBLIC Qt6::Core Qt6::Gui "-lmeadow" )
install( TARGETS qsimplex RUNTIME DESTINATION "bin" )
