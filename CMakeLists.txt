cmake_minimum_required(VERSION 3.14.0)
project(SmallPT VERSION 0.1.0)

include(verstring.cmake)
get_verstring(VERSTRING)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

set(BASE_FILES
    src/main.cpp
    src/smallpt.cpp
)

install(FILES ${CMAKE_BINARY_DIR}/Emu68.img DESTINATION .)

add_link_options(-static)
add_compile_options(-m68020 -m68881 -fomit-frame-pointer -fno-exceptions -ffreestanding -pedantic -pedantic-errors -Wall -Wextra -Werror)

add_executable(SmallPT
    ${BASE_FILES}
)

target_compile_definitions(SmallPT PRIVATE VERSION_STRING="${VERSTRING}")
