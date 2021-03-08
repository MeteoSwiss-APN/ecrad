#!/bin/bash
NVFORTRAN=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin/nvfortran
CUDA_HOME=/usr/local/cuda-10.2/ ${NVFORTRAN} -acc -mp=nonuma -O3 -gpu=cuda10.2 -Mcuda -Minfo=accel gpu_demo.F90 -o gpu_demo && ./gpu_demo

