#!/bin/bash

###############################
#Description
###############################
PROGRAM=Nevermore.sh
VERSION=0.2
AUTHOR="Minmin, Xue, PhD"
YEAR=2017
AFFILIATION="
University of Groningen, Nijenborgh 7, 9747AG Groningen, The Netherlands\n
Nanjing Universtiy of Aeronautics and Astronautics, 29 Yudao Str, Nanjing, China\n"

DESCRIPTION="\
# In order to save time and make things automatically run. I am trying to write this script.
# ############################ Minmin 28-2-2017
# This script is used to build a simulation box containing certain number of SMA copolymers and 
# a box of some kinds of lipids. Then put them together. 
# > 1. Building simulation systems
# > 2. Do minimization and equalibration in nvt ensemble with restraints
# > 3. Do equalibration in npt ensembles with different timesteps
# > 4. Production run"

echo "Program: \t $PROGRAM \n $DESCRIPTION \n Author:$AUTHOR Year: $YEAR \n $AFFILIATION \n"

##############################
# Gromacs version 5.x with command form "gmx command"
##############################
source /usr/local/gromacs-5.0.7/bin/GMXRC

GMXGMP="gmx grompp"
MDRUN="gmx mdrun "

EMMDP=~/common_files/em.mdp			# Energy minimization mdp file
EQ1NVT=~/common_files/eq1nvt.mdp		# pr-eq-NVT mdp file
EQ1NPT=~/common_files/eq1npt.mdp		# pr-eq-NPT mdp file
EQ2=~/common_files/eq2.mdp			# eq mdp file
MD=~/common_files/md.mdp			# md mdp file

RM=/bin/rm
CP=/bin/cp
MV=/bin/mv

$CP $EMMDP .
sed -i 's/TIMESTEP/0.001/g' ./em.mdp
sed -i 's/NSTEPS/1000/g' ./em.mdp

##############################
# Running parts
##############################

# EM ----------------------------------------------------------------
$GMXGMP -f em.mdp -o SYSTEM-EM.tpr -n SYSTEM.ndx -p SYSTEM-EM.top -maxwarn 10 -c SYSTEM.gro &> MD.log
$MDRUN  -deffnm SYSTEM-EM -v -stepout 1000 -maxh 72 -c SYSTEM-EMed.gro &>> MD.log

$CP $EQ1NVT .
sed -i 's/TIMESTEP/0.001/g' ./eq1nvt.mdp
sed -i 's/NSTEPS/200000/g' ./eq1nvt.mdp
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq1nvt.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq1nvt.mdp
sed -i 's/REFT/330 330 330/g' ./eq1nvt.mdp

$GMXGMP -f eq1nvt.mdp -o SYSTEM-EQ1NVT.tpr -n SYSTEM.ndx -p SYSTEM-POSRES.top -maxwarn 10 -c SYSTEM-EMed.gro &>> MD.log

$MDRUN  -deffnm SYSTEM-EQ1NVT -v -stepout 1000 -maxh 72 &>> MD.log

$CP $EQ1NPT .
sed -i 's/TIMESTEP/0.002/g;s/NSTEPS/500000/g' ./eq1npt.mdp
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq1npt.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq1npt.mdp
sed -i 's/REFT/330 330 330/g' ./eq1npt.mdp

$GMXGMP -f eq1npt.mdp -o SYSTEM-EQ1NPT.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ1NVT.cpt &>> MD.log

$MDRUN  -deffnm SYSTEM-EQ1NPT -v -stepout 1000 -maxh 72 &>> MD.log

$CP $EQ2 .
sed -i 's/TIMESTEP/0.005/g;s/NSTEPS/2000000/g' ./eq2.mdp
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq2.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq2.mdp
sed -i 's/REFT/330 330 330/g' ./eq2.mdp

$GMXGMP -f eq2.mdp -o SYSTEM-EQ2.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ1NPT.cpt &>> MD.log

$MDRUN  -deffnm SYSTEM-EQ2 -v -stepout 1000 -maxh 72 &>> MD.log

$CP $MD .
sed -i 's/TIMESTEP/0.015/g;s/NSTEPS/80000000/g' ./md.mdp
sed -i 's/GROUPS/DPPC SMA SOL/g' ./md.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./md.mdp
sed -i 's/REFT/330 330 330/g' ./md.mdp

$GMXGMP -f md.mdp -o SYSTEM-MD.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ2.cpt &>> MD.log

$MDRUN  -deffnm SYSTEM-MD -v -stepout 1000 -maxh 72 &>> MD.log

$RM \#* -f
