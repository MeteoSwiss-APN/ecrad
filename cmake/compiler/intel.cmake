if(__claw_base_compiler)
  return()
endif()
set(__claw_base_compiler YES)

# For Intel compiler
set(FC_BASIC_FLAGS "-convert big_endian -fpe0 -fp-model precise -ftz -fp-speculation safe")
set(FC_FPPFLAGS "") 
set(FC_RELEASE_OPTFLAGS "-O3")
set(FC_DEBUG_OPTFLAGS "-O0")
set(FC_RELEASE_WARNFLAGS "-warn all")
set(FC_DEBUG_WARNFLAGS "${FC_RELEASE_WARNFLAGS}")
set(FC_RELEASE_DEBUGFLAGS "-g")
set(FC_DEBUG_DEBUGFLAGS "${FC_RELEASE_DEBUGFLAGS}")
set(FC_OMPFLAG "-qopenmp -qopenmp-lib=compat")
set(FC_MODFLAG "-module")
set(CMAKE_EXE_LINKER_FLAGS_INIT "-lrt")

if(GPU_ACC)
    message( FATAL_ERROR "OpenACC not supported with Intel Compiler")
endif()


