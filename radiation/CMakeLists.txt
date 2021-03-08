file(GLOB radiation_SRC CONFIGURE_DEPENDS "*.F90")

target_sources(radiation PRIVATE ${radiation_SRC})

