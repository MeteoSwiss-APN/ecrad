if(__claw_base_compiler)
  return()
endif()
set(__claw_base_compiler YES)

# For PGI compiler
set(FC_BASIC_FLAGS "-Mbyteswapio -fpic")
set(FC_FPPFLAGS "-Mpreprocess") 
set(FC_RELEASE_OPTFLAGS "-O3")
set(FC_DEBUG_OPTFLAGS "-O0")
set(FC_RELEASE_WARNFLAGS "-Minform=inform")
set(FC_DEBUG_WARNFLAGS "${FC_RELEASE_WARNFLAGS}")
set(FC_RELEASE_DEBUGFLAGS "-g")
set(FC_DEBUG_DEBUGFLAGS "${FC_RELEASE_DEBUGFLAGS}")
set(FC_OMPFLAG "-mp=nonuma")
set(FC_MODFLAG "-J")

if("${GPU_ACC_TARGET_CC}" STREQUAL "")
set(GPU_TARGET_CC "")
else()
set(GPU_TARGET_CC ":${GPU_ACC_TARGET_CC}")
endif()

if(GPU_ACC)
    set(FC_BASIC_FLAGS "${FC_BASIC_FLAGS} -acc -ta=tesla${GPU_TARGET_CC},maxregcount:${GPU_ACC_TARGET_MAXREGCOUNT},rdc -Mcuda=ptxinfo -Minfo=accel -Mipa")
endif()

