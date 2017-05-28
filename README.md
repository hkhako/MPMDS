# MPMDS
Multi-parallel Molecular Dynamic Simulator

## Introduction
It is an atomistic molecular dynamic simulator which employs MPI (message passing interface) for massive parallel distributed computation.  The algorithm employs EAM (embedded atomic method) as a mean for force evaluation.  

## Build
The project was setup to compile with mpif90, but is tested with intel fortran compiler.  A compile script for mpif90 is included.
Just run ./compile.sh to compile  
