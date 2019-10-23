#!/bin/bash

#SBATCH -J POP2-Ca
#SBATCH -n 128
#SBATCH -t 72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=j.c.ramirez.palacios@student.rug.nl

# Set flags for mdrun (-cpi is set by the script when needed, no need to to include it here)
FLAGS='-dlb yes -maxh 72 -noappend -v'
# -cpo state.cpt  # needed if using -defmm

# Set minimum job run time (10 mn) - if job length less than $MINIMUMTIME don't continue next job
MINIMUMTIME=$(echo "10*60*1" | bc)

# Decrease maximum number of iterations ($maxrun arg passed to script vi SBATCH with --export=maxrun=$maxrun)
let maxrun--

# Set working directory
echo "The jobs are ran here:"
pwd

# Set active modules
# module add gcc/4.9.1 
# module add openmpi/gcc/1.8.1
# module add gromacs/openmpi/gcc/4.6.1
module add GROMACS/5.0.4-ictce-7.2.4-hybrid
# Start next job, but on hold waiting for this one
if [ $maxrun -gt 0 ] 
then
    sbatch --export=ALL,maxrun=$maxrun --dependency=afterok:$SLURM_JOB_ID launch_gromacs.sh >> job.id
    echo "Next job sbatch --export=ALL,maxrun=$maxrun --dependency=afterok:$SLURM_JOB_ID launch_gromacs.sh >> job.id"
fi

# Start mdrun job
STARTTIME=`date +%s`
# grompp -f mdrun.mdp -c equilibration.gro -p system.top -n index.ndx
if [ -f PREV.status ]
then
    STATUS_PREVIOUS=`tail -1 PREV.status`
    rm -f PREV.status
    if [ -f state.cpt -a -f state_prev.cpt ]
    then 
        cp state.cpt state_old_$STARTTIME.cpt
        cp state_prev.cpt state_prev_old_$STARTTIME.cpt
        if [ $STATUS_PREVIOUS == OK ]
        then
            echo "Run mdrun -cpi state.cpt"
            # srun mdrun_mpi -cpi state.cpt $FLAGS >> mdrun.log 2>&1
            srun mdrun_mpi -s topol.tpr -cpi state.cpt >> mdrun.log 2>&1

        elif [ $STATUS_PREVIOUS == NOTOK ] 
        then 
            echo "Run mdrun -cpi state_prev.cpt"
            srun mdrun_mpi -s topol.tpr -cpi state_prev.cpt >> mdrun.log 2>&1
        else
            echo "Status not recognized - input needed"
        fi
    else
        echo "ERROR: state.cpt and state_prev.cpt files not found"
    fi
elif [ ! -f state.cpt ]
then
    echo "Run mdrun"
    srun mdrun_mpi $FLAGS >> mdrun.log 2>&1
#    sleep 20
#    echo "after sleep run"
fi
MDRUNEXITCODE=$?

# If job terminated to early, kill next job
if [ $(echo "`date +%s`-$STARTTIME" | bc) -lt $MINIMUMTIME  -a  $maxrun -gt 0 ]
then
    echo "WARNING: job exited before min time, exit not ok"
    echo EARLY-STOP > PREV.status
    exit 1
elif [ $MDRUNEXITCODE != 0 ]
then
    echo NOTOK > PREV.status
else
    echo OK > PREV.status
fi
# test modify
