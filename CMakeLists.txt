cmake_minimum_required(VERSION 3.14.0)
project(SmallPT VERSION 0.1.1)

include(verstring.cmake)
get_verstring(VERSTRING)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

set(BASE_FILES
    src/main.c
)

install(FILES ${CMAKE_BINARY_DIR}/Emu68.img DESTINATION .)

add_link_options(-nostdlib -nostartfiles -static)
add_compile_options(-fno-exceptions -ffreestanding -pedantic -pedantic-errors -Wall -Wextra -Werror)

add_executable(SmallPT
    ${BASE_FILES}
)

target_compile_definitions(Emu68.elf PRIVATE VERSION_STRING="${VERSTRING}")
