#!/bin/bash
#SBATCH -J invspec            # job name
#SBATCH -o lslog.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested
#SBATCH -n 1               # total number of mpi tasks requested
#SBATCH -p normal      # queue (partition) -- normal, development, etc.
#SBATCH -t 03:00:00         # run time (hh:mm:ss) - 1.5 hours

#ssh login1 "mail -s 'SLURM Jobid=${SLURM_JOBID} Name=${SLURM_JOB_NAME} Began on Lonestar 5.' \
#	fl3537@utexas.edu <<< 'Job Started' "

module load matlab
#matlab -nodesktop -nosplash < $mname
matlab -nodesktop -nosplash < test_script.m

#ssh login1 "mail -s 'SLURM Jobid=${SLURM_JOB_ID} Name=${SLURM_JOB_NAME} Ended on Lonestar 5.' \
#	fl3537@utexas.edu  <<< 'Job Ended' "
