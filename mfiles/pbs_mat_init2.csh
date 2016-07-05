#!/bin/csh
#  Request maximum wallclock time for job
#PBS -l walltime=0:10:00
#  Total cores requested = <number of nodes> X <MPI procs/node>
#PBS -l select=1:ncpus=32:mpiprocs=1
# Specify how MPI processes are distributed on nodes
#PBS -l application=matlab
#  Request job name
#PBS -N init_pcui2
#  Request PBS job queue for job
##PBS -q background
##PBS -q standard
#PBS -q debug
#  Specify keep stdout/stderr files from job
#PBS -k eo
#  Indicate required  matlab licenses
#PBS -l MATLAB=1
#  Email
##PBS -m abe -M knelson3@stanford.edu 
# Indicate Project ID
#PBS -A ONRDC27755267
#  Request environment variables be exported from script
#PBS -V


## End of preamble, beginning of shell script ##
set JOBID=`echo $PBS_JOBID | cut -f1 -d.`
echo job $JOBID starting at `date` on `hostname`
cd $PBS_O_WORKDIR
echo starting in `pwd`

set outfile=init.out
set outdir=$PBS_O_WORKDIR
set indir=`pwd`

set infile=initialize_pcui.m

ls -l

module load matlab

echo matlab  $infile >& $outfile
matlab -nodesktop -nodisplay < $infile >& $outfile

set st=$status
echo execution ended at `date` with status $st

#chmod -R 755 $outdir
exit $st

