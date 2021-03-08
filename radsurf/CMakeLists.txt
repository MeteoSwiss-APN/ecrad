set(radsurf_SRC radsurf_properties.F90 radsurf_flux.F90 radsurf_intermediate.F90 radsurf_save.F90)

PREPEND(radsurf_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${radsurf_SRC})

target_sources(radsurf PRIVATE ${radsurf_SRC})

