#!/bin/csh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l application=matlab
#PBS -l MATLAB=1
#PBS -l MATLAB_Distrib_Comp_Engine=1 
#PBS -N Quad2
#PBS -q standard
#PBS -A ONRDC27755267
#PBS -l ccm=1

module load matlab/9.2.0
ccmrun matlab -nodisplay -nosplash -r "run('~/repository/pcui-3d/mfiles/QuadrantAnalysis2.m');" > ~/repository/pcui-3d/mfiles/quad2.out


