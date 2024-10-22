set(HISTORY_SOURCES
  hist_api.F90
  hist_buffer.F90
  hist_field.F90
)

add_library(history ${HISTORY_SOURCES})

target_compile_options(history PRIVATE -ffree-line-length-none)

add_subdirectory(util)
add_subdirectory(hash)

target_include_directories(history PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
