#!/bin/bash
NVFORTRAN=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin/nvfortran
#export PGI_ACC_DEBUG=1
#export PGI_ACC_SYNCHRONOUS=1
CUDA_HOME=/usr/local/cuda-10.2/ ${NVFORTRAN} -acc -Mcuda=rdc -mp=nonuma -Minfo -O3 -gpu=cuda10.2 mempool_demo.F90 -o mempool_demo_gpu && ./mempool_demo_gpu


