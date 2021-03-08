set(ecrad_SRC ecrad_driver_read_input.F90 ecrad_driver_config.F90 ecrad_driver.F90)

PREPEND(ecrad_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${ecrad_SRC})

set(ecrad_SRC ${ecrad_SRC} PARENT_SCOPE)
 


