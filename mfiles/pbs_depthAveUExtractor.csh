#!/bin/csh
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=32:mpiprocs=1
#PBS -l place=scatter
#PBS -l application=matlab
#PBS -N mean_CRe2_7
#PBS -q background
##PBS -q standard
##PBS -q debug
#PBS -k eo
#PBS -l MATLAB=1
#PBS -m abe -M knelson3@stanford.edu 
#PBS -A ONRDC27755267
#PBS -V


## End of preamble, beginning of shell script ##


cd $PBS_O_WORKDIR
set JOBID=`echo $PBS_JOBID | cut -f1 -d.`
echo job $JOBID starting at `date` on `hostname`
echo starting in `pwd`

set outfile=UmeanMat.out
set indir=`pwd`

set infile=~/mfiles/meanVelocityOnly.m


module load matlab

echo matlab  $infile >& $outfile
matlab -nodesktop -nodisplay < $infile >& $outfile

set st=$status
echo execution ended at `date` with status $st

exit $st
