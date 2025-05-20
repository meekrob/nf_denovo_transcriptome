#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --job-name=nf-start
#SBATCH --partition=week-long-cpu
#SBATCH --output=%x.%j.log # gives jobname.ID.log
# Available partitions
# day-long-cpu
# day-long-gpu
# day-long-highmem
# exp-gpu
# short-cpu*
# short-gpu
# short-highmem
# week-long-cpu
# week-long-gpu
# week-long-highmem
set -ue
set -o pipefail
datetime=$(date +"%Y-%m-%d_%H:%M:%S")
echo "[$datetime] $SLURM_JOB_NAME $@" # log the command line

mkdir -p $HOME/tmp
export TMP=$HOME/tmp
export TMPDIR=$TMP

cmd="nextflow run main.nf --profile slurm $@"
echo $cmd
time eval $cmd
