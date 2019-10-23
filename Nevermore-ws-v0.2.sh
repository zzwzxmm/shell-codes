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

GMXIRT="gmx insert-molecules"
GMXSOL="gmx solvate"
GMXION="gmx genion"
GMXNDX="gmx make_ndx"
GMXGMP="gmx grompp"
GMXEDF="gmx editconf"
MDRUN="gmx mdrun "

##############################
# Building parts -- You need a Directory contain all itp files and executable files
##############################
DEPENDENCIES=(insane.py)

# ----------------------------
# Top file headers
# Martini itp files, version 2.2

MARTINI="#include \"/home/minmin/common_files/martini_v2.2p.itp\" "
MARTINIEM="#include \"/home/minmin/common_files/martini_v2.2pem.itp\" "
MARTINIPARTICLES="#include \"/home/minmin/common_files/martini_particles.itp\" "

# Martini lipid itp files
DPPC="#include \"/home/minmin/common_files/martini_DPPC.itp\" "
DOPC="#include \"/home/minmin/common_files/martini_DOPC.itp\" "
DLPC="#include \"/home/minmin/common_files/martini_DLPC.itp\" "
DTPC="#include \"/home/minmin/common_files/martini_DTPC.itp\" "

# SMA itp file
SMA="#include \"/home/minmin/common_files/sma_cg-23units.itp\" "

# Position Restraints files
SMAPOSRES="#include \"/home/minmin/common_files/sma-posres.itp\" "	# If SMAPOSRES is used, SMA should be abandoned
DPPCPOSRES="\
#ifdef POSRES
#include \"/home/minmin/common_files/dppc-posres.itp\"
#endif
"
DOPCPOSRES="\
#ifdef POSRES
#include \"/home/minmin/common_files/dopc-posres.itp\"
#endif
"
DLPCPOSRES="\
#ifdef POSRES
#include \"/home/minmin/common_files/dlpc-posres.itp\"
#endif
"
DTPCPOSRES="\
#ifdef POSRES
#include \"/home/minmin/common_files/dtpc-posres.itp\"
#endif
"

# -----------------------------------
# Starting step, variables

# ----------------------------------
# System description
#               -----------------
#               |     PW box    |
#               -----------------
#               |    SMA box    |
#               -----------------
#               |   Lipid box   |
#               -----------------
#-----------------------------------

# Variables ------------------------------------------------------------
LIPIDS=(DPPC)					# e.g. (DPPC:1 DOPC:2 DLPC:3);
APL=0.65					# Area per lipid, depend on lipid types. 
SALTCON=0.15					# salt concentration;
LBOX=(20 20 10)					# lipids boxsize in nm;
BOXTYPE="cubic"					# boxtype;
SOLVENT="PW"					# Solvent type, W or PW 
SBOX=(20 20 5)					# SMA box size
SNUM=10						# SMA copolymer number
PBOX=(20 20 5)					# PW box size

EMMDP=~/common_files/em.mdp			# Energy minimization mdp file
EQ1NVT=~/common_files/eq1nvt.mdp		# pr-eq-NVT mdp file
EQ1NPT=~/common_files/eq1npt.mdp		# pr-eq-NPT mdp file
EQ2=~/common_files/eq2.mdp			# eq mdp file
MD=~/common_files/md.mdp			# md mdp file
SMAGRO=~/common_files/sma_cg-23units.gro	# SMA gro structure file
INSANE=~/common_files/insane.py			# Insane.py to build lipids
PWATER=~/common_files/pwater.gro		# Polarizable water structure gro file

RM=/bin/rm
CP=/bin/cp
MV=/bin/mv

$CP $EMMDP .
sed -i 's/TIMESTEP/0.001/g' ./em.mdp
sed -i 's/NSTEPS/1000/g' ./em.mdp

# Lipids -----------------------------------------------------------------
LipidOption=$(for ((i=0; i<${#LIPIDS[@]}; i++)); do echo -n "-l ${LIPIDS[$i]} "; done)
$INSANE $LipidOption -salt $SALTCON -x ${LBOX[0]} -y ${LBOX[1]} -z ${LBOX[2]} -d 0 -a $APL -pbc ${BOXTYPE} -sol ${SOLVENT} -o LIPIDS.gro -p LIPIDS.top
# insane builds the box from 0 0 0 to 20 20 10; Zeta direction thickness depends on lipid type.

# create sma box ---------------------------------------------------------
$GMXIRT -ci $SMAGRO -o SMA.gro -try 500 -nmol $SNUM -box ${SBOX[@]} &> 1-GMXIRT-SMA.log
# SMA box. g_genbox before 5.0. after 5.0, genbox is splited into solvate and insert-molecules.

# write sma top file -----------------------------------------------------
read -a SMAarray <<< $(echo $(cat 1-GMXIRT-SMA.log | grep 'Output configuration contains'))
SmaNumber=${SMAarray[6]}
sh -c "cat > SMA.top" << EOF
$MARTINI
$SMA
$MARTINIPARTILES

[ system ]
; name
SMA-SYSTEM

[ molecules ]
; name  number
SMA-Copolymer	  $((SmaNumber/23))
EOF

# solvate sma gro file ---------------------------------------------------
$GMXSOL -cp SMA.gro -cs $PWATER -o SMA-PW.gro -p SMA.top &> 1-GMXSOL-1.log
# write top file for sma-water box ---------------------------------------
read -a PWarray <<< $(echo $(cat 1-GMXSOL-1.log | grep 'Generated solvent containing'))
PwNumber=${PWarray[6]}
sh -c "cat >> SMA.top" << EOF
PW        $((PwNumber))
EOF

# make ndx file of sma-pw ------------------------------------------------
echo -en "q\n" | $GMXNDX -f SMA-PW.gro -o SMA-PW.ndx  &> 1-GMXNDX-SMA-PW.log

# make sma-pw-em.tpr file for genion -------------------------------------
$GMXGMP -f em.mdp -o SMA-PW-EM.tpr -n SMA-PW.ndx -p SMA.top -maxwarn 10 -c SMA-PW.gro &> 1-GMXGMP-SMA-PW.log

# genion in sma-pw-ion.gro file ------------------------------------------
echo -en "7\n" | $GMXION -s SMA-PW-EM.tpr -n SMA-PW.ndx -o SMA-PW-ION.gro -p SMA.top -conc $SALTCON -neutral  &> 1-GMXION-SMA-PW.log

# create a water box up the sma ------------------------------------------
$GMXSOL -cs $PWATER -o PW.gro -box ${PBOX[@]} &> 1-GMXSOL-PW.log

# write pw top file ------------------------------------------------------
read -a PWarray <<< $(echo $(cat 1-GMXSOL-PW.log | grep 'Output configuration contains'))
PwNumber=${PWarray[6]}
sh -c "cat > PW.top" << EOF
$MARTINI
$MARTINIPARTICLES

[ system ]
; name
PW-SYSTEM

[ molecules ]
; name  number
PW     $PwNumber
EOF

# make ndx file for pw.gro -----------------------------------------------
echo -en "q\n" | $GMXNDX -f PW.gro -o PW.ndx &> 1-GMXNDX-PW.log

# create pw-em.tpr for gmx genion ----------------------------------------
$GMXGMP -f em.mdp -o PW-EM.tpr -n PW.ndx -p PW.top -maxwarn 10 -c PW.gro &> 1-GMXMDP-PW-EM.log

# create pw-ion gro file -------------------------------------------------
echo -en "2\n" | $GMXION -s PW-EM.tpr -n PW.ndx -o PW-ION.gro -p PW.top -conc $SALTCON -neutral  &> 1-GMXION-PW.log

# merge the three files in one -------------------------------------------
$GMXEDF -f LIPIDS.gro -o LIPIDS.pdb &> 1-GMXEDF-LIPID.log
$GMXEDF -f SMA-PW-ION.gro -o SMA-PW-ION.pdb -translate 0 0 ${LBOX[2]} &> 1-GMXEDF-SMA-PW-ION.log
$GMXEDF -f PW-ION.gro -o PW-ION.pdb -translate 0 0 $((LBOX[2]+SBOX[2])) &> 1-GMXEDF-PW-ION.log
cat LIPIDS.pdb SMA-PW-ION.pdb PW-ION.pdb | grep "^ATOM" > SYSTEM.pdb
$GMXEDF -f SYSTEM.pdb -o SYSTEM.gro -box ${LBOX[0]} ${LBOX[1]} $((LBOX[2]+SBOX[2]+PBOX[2])) -angles 90 90 90 &> 1-GMXEDF-SYSTEM.log

# create top file for the system -----------------------------------------
cat LIPIDS.top SMA.top PW.top > TEMP.top
for i in { "#" "\[" ";" "^$" "INSNANE" "SYSTEM" }; do sed -i "/$i/d" TEMP.top; done
sh -c "cat > SYSTEM.top" << EOF
$MARTINI
$DPPC
$SMA
$MARTINIPARTICLES

[ system ]
; name
SYSTEM

[ molecules ]
; name    number
`cat TEMP.top`
EOF

sh -c "cat > SYSTEM-EM.top" << EOF
$MARTINIEM
$DPPC
$SMA
$MARTINIPARTICLES

[ system ]
; name
SYSTEM

[ molecules ]
; name    number
`cat TEMP.top`
EOF

sh -c "cat > SYSTEM-POSRES.top" << EOF
$MARTINI
$DPPC
$DPPCPOSRES
$SMAPOSRES
$MARTINIPARTICLES

[ system ]
; name
SYSTEM

[ molecules ]
; name    number
`cat TEMP.top`
EOF

# creater ndx file for the system ----------------------------------------
echo -en "6|7|8|9|10\nname 25 SMA\n4|11\nname 26 SOD\n5|12\nname 27 CLA\n3|26|27\nname 28 SOL\nq\n" | $GMXNDX -f SYSTEM.gro -o SYSTEM.ndx
#or use expect ( not recommended ) or gmx select

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
