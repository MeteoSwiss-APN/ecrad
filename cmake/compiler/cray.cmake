if(__claw_base_compiler)
  return()
endif()
set(__claw_base_compiler YES)

# For Cray compiler
set(FC_BASIC_FLAGS "-emf -fpic -hbyteswapio -hfp_trap")
set(FC_FPPFLAGS "-eT") 
set(FC_RELEASE_OPTFLAGS "-hflex_mp=conservative -hfp1 -hadd_paren -hnocaf -hvector3 -O3")
set(FC_DEBUG_OPTFLAGS "-O0 -hfp1")
set(FC_RELEASE_WARNFLAGS "")
set(FC_DEBUG_WARNFLAGS "${FC_RELEASE_WARNFLAGS}")
set(FC_RELEASE_DEBUGFLAGS "")
set(FC_DEBUG_DEBUGFLAGS "${FC_RELEASE_DEBUGFLAGS}")
set(FC_OMPFLAG "-homp")
set(FC_MODFLAG "-J")
set(CMAKE_EXE_LINKER_FLAGS_INIT "-Ktrap=fp")

if(GPU_ACC)
    message( FATAL_ERROR "OpenACC not supported with Cray Compiler")
endif()


