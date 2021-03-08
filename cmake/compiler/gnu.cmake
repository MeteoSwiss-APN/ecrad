if(__claw_base_compiler)
  return()
endif()
set(__claw_base_compiler YES)

# For GNU compiler
set(FC_BASIC_FLAGS "-fconvert=big-endian -fno-range-check")
set(FC_FPPFLAGS "-cpp") 
set(FC_RELEASE_OPTFLAGS "-O3")
set(FC_DEBUG_OPTFLAGS "-O0")
set(FC_RELEASE_WARNFLAGS "-Wall -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable")
set(FC_DEBUG_WARNFLAGS "-Wall")
set(FC_RELEASE_DEBUGFLAGS "-g -ffpe-trap=invalid,zero,overflow")
set(FC_DEBUG_DEBUGFLAGS "-g -ffpe-trap=invalid,zero,overflow -fcheck=bounds -finit-real=snan -fbacktrace")
set(FC_OMPFLAG "-fopenmp")
set(FC_MODFLAG "-J")

if(GPU_ACC)
    message( FATAL_ERROR "OpenACC not supported with GNU Compiler")
endif()

