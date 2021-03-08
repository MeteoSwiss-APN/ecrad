#!/bin/bash
NVFORTRAN=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin/nvfortran
${NVFORTRAN} -mp=nonuma -O3 mempool_demo.F90 -o mempool_demo && ./mempool_demo

