#!/bin/sh
#PBS -N channel1p-equil-E0.002
#PBS -e channel1p-equil-E0.002.err
#PBS -o channel1p-equil-E0.002.out
#PBS -l nodes=6:ppn=16
#PBS -l walltime=120:00:00

echo "Working dir: " $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo "Nodes used: " `uniq $PBS_NODEFILE`
NPROCS=`wc -l < $PBS_NODEFILE`
cp $PBS_NODEFILE test.nodefile
echo "Cores used: " $NPROCS
echo "Starting at: " `date`

/opt/openmpi/bin/mpirun -machinefile $PBS_NODEFILE -np $NPROCS ./sedi

echo "Ending at: " `date`
