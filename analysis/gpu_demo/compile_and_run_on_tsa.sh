#!/bin/bash
pgfortran -acc -mp=nonuma -O3 -g -ta=tesla:cc70,maxregcount:32,rdc -Mcuda=ptxinfo -Minfo=accel gpu_demo.F90 -o gpu_demo && ./gpu_demo

