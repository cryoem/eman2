#!/bin/bash
set echo on

#$ -pe mpi XXX_SXMPI_NPROC_XXX
#$ -N XXX_SXMPI_JOB_NAME_XXX
#$ -cwd 

# Environment
source ~/.bash_profile

#Hostfiles
set hosts = `cat $PE_HOSTFILE | cut -f 1 -d . | sort | fmt -w 30 | sed 's/ /,/g'`
#echo $hosts > $dir/hosts.$JOB_ID
cat $PE_HOSTFILE | awk '{print $1}'  > hostfile.$JOB_ID

mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX