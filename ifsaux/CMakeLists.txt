file(GLOB ifsaux_SRC CONFIGURE_DEPENDS "*.F90")

target_sources(ifsaux PRIVATE ${ifsaux_SRC})
