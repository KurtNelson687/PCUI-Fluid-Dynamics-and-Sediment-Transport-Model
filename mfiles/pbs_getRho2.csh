#!/bin/csh
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=32:mpiprocs=1
#PBS -l application=matlab
#PBS -N getRho200_11cup
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

set outfile=getRho02.out
set indir=`pwd`

set infile=~/repository/pcui-3d/mfiles/getRhoPrimeMean2.m


module load matlab

echo matlab  $infile >& $outfile
matlab -nodesktop -nodisplay < $infile >& $outfile

set st=$status
echo execution ended at `date` with status $st

exit $st

