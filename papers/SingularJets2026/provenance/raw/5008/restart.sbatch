#!/bin/bash
#SBATCH --job-name=sbb5008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=06:00:00
#SBATCH --partition=genoa
#SBATCH --output=slurm-restart-%x-%j.out
#SBATCH --error=slurm-restart-%x-%j.err
cd "/gpfs/work2/0/nctt0620/vatsal/2026-07-03-Singular-Bursting-Bubbles-Bo0-OhSweep-L14/simulationCases/5008"
echo "RESTART 5008 from snapshot-0.490625  gas=3 keStop=1e6 threads=64  $(date)"
OMP_NUM_THREADS=64 OMP_PROC_BIND=spread OMP_PLACES=cores ./drill case.params
echo "FINISHED 5008 rc=$? $(date)"
