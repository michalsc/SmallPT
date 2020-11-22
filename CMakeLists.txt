cmake_minimum_required(VERSION 3.14.0)
project(SmallPT VERSION 0.1.0)

include(verstring.cmake)
get_verstring(VERSTRING)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

set(SUPPORTED_TARGETS "AmigaOS" "Qt")
set(TARGET "AmigaOS" CACHE STRING "One of target machines: ${SUPPORTED_TARGETS}")
set_property(CACHE TARGET PROPERTY STRINGS ${SUPPORTED_TARGETS})

set(TARGET_FILES "")
set(BASE_FILES
    src/smallpt.cpp
)

if (${TARGET} IN_LIST SUPPORTED_TARGETS)
    message("-- Selected target machine: ${TARGET}")
    if(${TARGET} STREQUAL "AmigaOS")
        list(APPEND TARGET_FILES src/main_amiga.cpp)
    elseif(${TARGET} STREQUAL "Qt")
        find_package(Qt5 COMPONENTS Widgets LinguistTools REQUIRED)
        target_link_libraries(SmallPT Qt5::Widgets)
        list(APPEND TARGET_FILES src/main_qt.cpp)
    endif()
else()
    message(FATAL_ERROR "Wrong target machine specified: ${TARGET}")
endif()

install(FILES ${CMAKE_BINARY_DIR}/SmallPT DESTINATION .)

add_compile_options(-m68020 -m68881 -fomit-frame-pointer -fno-exceptions -pedantic -pedantic-errors -Wall -Wextra -Werror)

add_executable(SmallPT
    ${TARGET_FILES}
    ${BASE_FILES}
)

target_compile_definitions(SmallPT PRIVATE VERSION_STRING="${VERSTRING}")
