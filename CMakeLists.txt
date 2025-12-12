cmake_minimum_required(VERSION 3.20)
# as your project grows you will need to add more source filenames and 
# library names in the project config section for the project to build
# wxWidgets settings should only need changing if you switch from static 
# to dynamic built

### project config
      # your project's name
      set(projectName "GrunObject")
      # names of source files in your project here (.cpp only)
      set(SOURCES 
		GrunObject.cpp
		GrunItem.cpp
      )
      # names of libraries your project needs to LINK to here (if ld.exe complains it can't find things, chances are you need to add a library name here from your lib/ folder)
      # SQLiteCpp depends on sqlite3!
      # C++23 std::println under Windows needs stdc++exp here
      set(projectLibraries
            #SQLiteCpp
            #sqlite3
            stdc++exp
      )
      set(localCPPDevRoot "C:/dev/cpp")
### end project config

### wxWidgets config - have you set Environment variable WXWIN? If not, do it now.
      # base directory where wxWidgets was installed to (I install wx to my "cpp dev root" which is C:/cpp/ and I set ENV{WXWIN} to that path)
      set(wxWidgets_ROOT_DIR $ENV{WXWIN})
      # Dynamic (OFF) or Static (ON) build?
      set(wxWidgets_USE_STATIC ON CACHE BOOL "Use static libraries" FORCE)
      # local machine/local build settings (look in your dev root's lib/ folder for folder names, use *_lib for static, *_dll for dynamic)
      set(localMachine_WX_LIB_NAME "gcc_x64_lib")
      # local machine specific configuration - the configuration name (look in lib/gcc_*/ folder for folder names)
      set(wxWidgets_CONFIGURATION "mswu")
      # the wx components the project needs !! NOTE: do not enclose in quotes/double quotes !!
      set(wxPackages gl core base)
### end wxWidgets config


######### you shouldn't have to configure anything else below here
project(${projectName} VERSION 0.2.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)

set(wxWidgets_LIB_DIR ${wxWidgets_ROOT_DIR}/lib/${localMachine_WX_LIB_NAME})
set(wxWidgets_INCLUDE_DIRS
      "${wxWidgets_ROOT_DIR}/include"
      "${wxWidgets_ROOT_DIR}/lib/${localMachine_WX_LIB_NAME}/${wxWidgets_CONFIGURATION}"
   )

include_directories(
      ${workspaceFolder}
      ${localCPPDevRoot}/include
      $ENV{WXWIN}/lib/${localMachine_WX_LIB_NAME}/${wxWidgets_CONFIGURATION}
      )

link_directories($ENV{WXWIN}/lib/${localMachine_WX_LIB_NAME} ${localCPPDevRoot}/lib)

# Find wxWidgets for static linking
find_package(wxWidgets 3.2 COMPONENTS ${wxPackages} REQUIRED)

include(${wxWidgets_USE_FILE} ${localCPPDevRoot}/include)

# these instructions are to build the library code in the project folder
# If there are source files, build a static library. If not, create a header-only INTERFACE library.
if(SOURCES)
    add_library(${projectName} STATIC ${SOURCES})
    target_include_directories(${projectName} PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${localCPPDevRoot}/include
        ${wxWidgets_INCLUDE_DIRS}
    )
    target_link_libraries(${projectName} PUBLIC ${wxWidgets_LIBRARIES} ${projectLibraries})
    target_compile_definitions(${projectName} PRIVATE ${wxWidgets_DEFINITIONS})
else()
    # Header-only / interface library (no object code produced)
    add_library(${projectName} INTERFACE)
    target_include_directories(${projectName} INTERFACE
        ${CMAKE_SOURCE_DIR}
        ${localCPPDevRoot}/include
        ${wxWidgets_INCLUDE_DIRS}
    )
    target_link_libraries(${projectName} INTERFACE ${wxWidgets_LIBRARIES} ${projectLibraries})
    target_compile_definitions(${projectName} INTERFACE ${wxWidgets_DEFINITIONS})
endif()

# Optional small test program to verify the library; only added if test/main.cpp exists.
option(BUILD_TEST_APP "Build small test executable to validate the library" ON)
if(BUILD_TEST_APP)
    set(TEST_MAIN "${CMAKE_SOURCE_DIR}/test/main.cpp")
    if(EXISTS ${TEST_MAIN})
        add_executable(${projectName}Test ${TEST_MAIN})
        target_link_libraries(${projectName}Test PRIVATE ${projectName})
        target_include_directories(${projectName}Test PRIVATE ${CMAKE_SOURCE_DIR})
    else()
        message(STATUS "No ${TEST_MAIN} found; skipping test executable creation.")
    endif()
endif()
# end instructions for bulding libraries in project folder

# commented out for building libraries only - add sources to executable
#add_executable(${projectName} WIN32 ${SOURCES})

# target_include_directories(${projectName} PRIVATE ${wxWidgets_INCLUDE_DIRS})
# target_link_libraries(${projectName} ${wxWidgets_LIBRARIES} ${projectLibraries})
# target_compile_definitions(${projectName} PRIVATE ${wxWidgets_DEFINITIONS})

# Handle project resources...
# Define the source and destination resource directories
# set(RESOURCE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/resources")
# set(RESOURCE_DEST_DIR "${CMAKE_BINARY_DIR}/resources")
# # Create a custom target to copy the resources recursively
# add_custom_target(copy_resources ALL
#     COMMAND ${CMAKE_COMMAND} -E copy_directory ${RESOURCE_SOURCE_DIR} ${RESOURCE_DEST_DIR}
#     COMMENT "Copying resource files recursively to the build directory..."
# )

# # Create a custom command to clean the destination directory
# add_custom_command(TARGET copy_resources PRE_BUILD
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${RESOURCE_DEST_DIR}
#     COMMENT "Removing outdated resource files from the build directory..."
# )

# # Add a custom command to copy the resources recursively
# add_custom_command(TARGET copy_resources POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_directory ${RESOURCE_SOURCE_DIR} ${RESOURCE_DEST_DIR}
#     COMMENT "Copying resource files recursively to the build directory..."
# )

# # Ensure the resources are copied before building the executable
# add_dependencies(${projectName} copy_resources)