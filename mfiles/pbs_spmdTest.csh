#!/bin/csh
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l application=matlab
#PBS -l MATLAB=1
#PBS -l MATLAB_Distrib_Comp_Engine=1 
#PBS -N spmdTest
#PBS -q standard
#PBS -A ONRDC27755267
#PBS -l ccm=1

set outfile=spmdTest.out
set infile=~/repository/pcui-3d/mfiles/spmdTest.m
module load matlab/9.2.0

#matlab -nodisplay -nosplash -r < $infile >& $outfile
ccmrun matlab -nodisplay -nosplash -r "run('~/repository/pcui-3d/mfiles/QuadrantAnalysis2.m');" > ~/repository/pcui-3d/mfiles/Quad2.out


