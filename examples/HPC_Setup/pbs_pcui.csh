#!/bin/csh

#PBS -N C3_1pAd2_bug
#PBS -q debug
#PBS -l walltime=00:10:00
##PBS -q standard
##PBS -l walltime=30:00:00
##PBS -q standard-long  
##PBS -l walltime=200:00:00
#PBS -V
#PBS -l select=15:ncpus=32:mpiprocs=32
#PBS -j oe
#PBS -r n
##PBS -m be -M knelson3@stanford.edu 
#PBS -A ONRDC27755267
#PBS -l place=scatter:excl

# PBS to LSF environmental variable equivalence settings

set input=cavity
set outfile=output.out

setenv JOBID `echo ${PBS_JOBID} | cut -f1 -d.`
echo "JOBID,PBS_JOBID=$JOBID,$PBS_JOBID"

setenv LOGNAME ${PBS_O_LOGNAME}
echo "LOGNAME,PBS_O_LOGNAME=$LOGNAME,$PBS_O_LOGNAME"

setenv NSLOTS "`wc -l ${PBS_NODEFILE} | cut -f1 -d' '`"
echo "**** NSLOTS=${NSLOTS}"

echo MPI job $JOBID starting at `date` on `hostname`
echo starting in `pwd`

# exports or setenv

setenv JOBID `echo $PBS_JOBID  | awk -F'.' '{print $1}'`

#setenv MPI_DSM_DISTRIBUTE

# For OpenMP jobs only! (NOT OPENMPI)
#setenv OMP_NUM_THREADS 1

## Get the execution host name
set exechost=`hostname -s`
echo exechost is $exechost

set jobhost=`uname -n`

# Copy Files to $TMPD which is usually /usr/var/tmp/userid/job #
cd $PBS_O_WORKDIR
cp $HOME/PCUI/sourceFiles/Domain3p2/* ./

make clean
make

ls -ltr
unlimit

#  The below section of commands verifies access to all the compute nodes to be used.
#========================================================================
setenv MY_HOSTS `cat ${PBS_NODEFILE} | tr ' ' '\n' `
setenv xTMP_HOSTS `echo ${MY_HOSTS} |  sed 's/ /\n/g' | sort -u `

foreach host ($xTMP_HOSTS)
  echo "Working on $host ...."
  /usr/bin/ssh -o StrictHostKeyChecking=no $host pwd
end
#========================================================================

#module unload mpi/sgi_mpi/2.03
#module load mpi/intelmpi-4.0.3
#module load compiler/intel12.1 mpi/intelmpi-4.0.3 

#intelmpirun.pbs ./${input} ${param_file} > ${outfile}
#mpirun -np 432 ./${input} > ${outfile}
aprun -n 480 ./${input} > ${outfile}

cp ${outfile} $PBS_O_WORKDIR

set st=$status
echo "MPI program ended with status $st on `date`"
exit $st

