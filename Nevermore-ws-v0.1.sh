#!/bin/bash
# Gromacs version 5.0.7
source /usr/local/gromacs-5.0.7/bin/GMXRC
MDRUN="gmx mdrun "

#########################################################################

PROGRAM=SMA_BUILDER.sh
VERSION=0.1
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
# > 1. Build two boxes of SMAs and lipids separately
# > 2. Put them together
# > 3. Do minimization and equalibration in nvt ensemble with restraints
# > 4. Do equalibration in npt ensembles with different timesteps
# > 5. Production run"

##############################
# Building parts -- You need a Directory contain all itp files and executable files
##############################
DEPENDENCIES=(insane.py)

MARTINI="\
#include \"/home/minmin/common_files/martini_v2.2p.itp\"
"
MARTINIEM="\
#include \"/home/minmin/common_files/martini_v2.2pem.itp\"
"
MARTINIPARTICLES="\
#include \"/home/minmin/common_files/martini_particles.itp\"
"
DPPC="\
#include \"/home/minmin/common_files/martini_DPPC.itp\"
"
DOPC="\
#include \"/home/minmin/common_files/martini_DOPC.itp\"
"
DLPC="\
#include \"/home/minmin/common_files/martini_DLPC.itp\"
"
DTPC="\
#include \"/home/minmin/common_files/martini_DTPC.itp\"
"
SMA="\
#include \"/home/minmin/common_files/sma_cg-23units.itp\"
"
SMAPOSRES="\
#include \"/home/minmin/common_files/sma-posres.itp\"
"
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

# Starting step 

# options:
LIPIDS=(DPPC)				# e.g. (DPPC:1 DOPC:2 DLPC:3);
SALTCON=0.15				# salt concentration;
LBOX=(20 20 10)				# lipids boxsize in nm;
BOXTYPE="cubic"				# boxtype;
SOLVENT="PW"
SBOX=(20 20 5)
SNUM=10
PBOX=(20 20 5)
EMMDP=~/common_files/em.mdp
EQ1NVT=~/common_files/eq1nvt.mdp
EQ1NPT=~/common_files/eq1npt.mdp
EQ2=~/common_files/eq2.mdp
MD=~/common_files/md.mdp
SMAGRO=~/common_files/sma_cg-23units.gro
INSANE=~/common_files/insane.py
PWATER=~/common_files/pwater.gro


GMXIRT="gmx insert-molecules"
GMXSOL="gmx solvate"
GMXION="gmx genion"
GMXNDX="gmx make_ndx"
GMXGMP="gmx grompp"

# Lipids -----------------------------------------------------------------
LipidOption=$(for ((i=0; i<${#LIPIDS[@]}; i++)); do echo -n "-l ${LIPIDS[$i]} "; done)
$INSANE $LipidOption -salt $SALTCON -x ${LBOX[0]} -y ${LBOX[1]} -z ${LBOX[2]} -d 0 -pbc ${BOXTYPE} -sol ${SOLVENT} -o LIPIDS.gro -p LIPIDS.top
# insane builds the box from 0 0 0 to 20 20 10; Zeta direction thickness depends on lipid type.
echo lipid box build
# ------------------------------------------------------------------------

# create sma box ---------------------------------------------------------
gmx insert-molecules -ci $SMAGRO -o SMA.gro -try 500 -nmol $SNUM -box ${SBOX[@]} &> INSERT.log
# SMA box. g_genbox before 5.0. after 5.0, genbox is splited into solvate and insert-molecules.
echo sma box build
# ------------------------------------------------------------------------

# write sma top file -----------------------------------------------------
read -a SMAarray <<< $(echo $(cat INSERT.log | grep 'Output configuration contains'))
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
echo sma toplogy write
# ------------------------------------------------------------------------

# solvate sma gro file ---------------------------------------------------
gmx solvate -cp SMA.gro -cs $PWATER -o SMA-PW.gro -p SMA.top &> SMA-PW.log
# write top file for sma-water box
read -a PWarray <<< $(echo $(cat SMA-PW.log | grep 'Generated solvent containing'))
PwNumber=${PWarray[6]}
sh -c "cat >> SMA.top" << EOF
PW        $((PwNumber))
EOF
echo sma-water box build
# ------------------------------------------------------------------------

# make ndx file of sma-pw ------------------------------------------------
gmx make_ndx -f SMA-PW.gro -o SMA-PW.ndx << EOF &> SMA-PW-NDX.log
q\r
EOF
#or use expect
echo sma-pw ndx file write
# ------------------------------------------------------------------------

# make sma-pw-em.tpr file for genion -------------------------------------
gmx grompp -f $EMMDP -o SMA-PW-EM.tpr -n SMA-PW.ndx -p SMA.top -maxwarn 10 -c SMA-PW.gro &> SMA-PW-MDP.log
# creater tpr file for gmx genion
# ------------------------------------------------------------------------

# genion in sma-pw-ion.gro file ------------------------------------------
gmx genion -s SMA-PW-EM.tpr -n SMA-PW.ndx -o SMA-PW-ION.gro -p SMA.top -conc $SALTCON -neutral << EOF &> SMA-PW-ION.log
7\r
EOF
echo sma-pw-ion gro file write
# ------------------------------------------------------------------------

# create a water box up the sma ------------------------------------------
gmx solvate -cs $PWATER -o PW.gro -box ${PBOX[@]} &> PW.log
echo pw gro file write
# ------------------------------------------------------------------------

# write pw top file ------------------------------------------------------
read -a PWarray <<< $(echo $(cat PW.log | grep 'Output configuration contains'))
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
echo pw top file write
# ------------------------------------------------------------------------

# make ndx file for pw.gro -----------------------------------------------
gmx make_ndx -f PW.gro -o PW.ndx << EOF &> SMA-PW-NDX.log
q\r
EOF
echo pw ndx file write
# ------------------------------------------------------------------------

# create pw-em.tpr for gmx genion ----------------------------------------
gmx grompp -f $EMMDP -o PW-EM.tpr -n PW.ndx -p PW.top -maxwarn 10 -c PW.gro &> PW-MDP.log
echo pw-em tpr file write
# ------------------------------------------------------------------------

# create pw-ion gro file -------------------------------------------------
gmx genion -s PW-EM.tpr -n PW.ndx -o PW-ION.gro -p PW.top -conc $SALTCON -neutral << EOF &> PW-ION.log
2\r
EOF
echo pw-ion gro file write
# ------------------------------------------------------------------------

# merge the three files in one -------------------------------------------
gmx editconf -f LIPIDS.gro -o LIPIDS.pdb &> EDITCONF.log
# 1. -center 0 0 0; 2. -translate
gmx editconf -f SMA-PW-ION.gro -o SMA-PW-ION.pdb -translate 0 0 ${LBOX[2]} &>> EDITCONF.log
# move up along zeta axis
gmx editconf -f PW-ION.gro -o PW-ION.pdb -translate 0 0 $((LBOX[2]+SBOX[2])) &>> EDITCONF.log
# move up along zeta axis

cat LIPIDS.pdb SMA-PW-ION.pdb PW-ION.pdb | grep "^ATOM" > SYSTEM.pdb
editconf -f SYSTEM.pdb -o SYSTEM.gro -box ${LBOX[0]} ${LBOX[1]} $((LBOX+SBOX+PBOX)) -angles 90 90 90 &>> EDITCONF.log
echo simulation system build
# ------------------------------------------------------------------------

# create top file for the system -----------------------------------------
cat LIPIDS.top SMA.top PW.top > TEMP.top
sed -i '/#/d' TEMP.top
sed -i '/;/d' TEMP.top
sed -i '/\[/d' TEMP.top
sed -i '/^$/d' TEMP.top
sed -i '/INSANE/d' TEMP.top
sed -i '/SYSTEM/d' TEMP.top
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
echo system top file write

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
echo system top file write

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
echo system-posres top file write
# ------------------------------------------------------------------------

# creater ndx file for the system ----------------------------------------
expect -c "
spawn gmx make_ndx -f SYSTEM.gro -o SYSTEM.ndx
expect {
\"> \" {
	send -- \"6|7|8|9|10\r\"
	send -- \"name 25 SMA\r\"
	send -- \"4|11\r\"
	send -- \"name 26 SOD\r\"
	send -- \"5|12\r\"
	send -- \"name 27 CLA\r\"
	send -- \"3|26|27\r\"
	send -- \"name 28 SOL\r\"
	send -- \"q\r\"
	}
}
interact
"
#or use expect
echo system ndx file write
# ------------------------------------------------------------------------

# create tpr file for the system EM --------------------------------------
gmx grompp -f $EMMDP -o SYSTEM-EM.tpr -n SYSTEM.ndx -p SYSTEM-EM.top -maxwarn 10 -c SYSTEM.gro 
echo system tpr em file write 
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------

##############################
# Running parts
##############################

# Launch the parallel job
echo $MDRUN
$MDRUN  -deffnm SYSTEM-EM -v -stepout 1000 -maxh 72 -c SYSTEM-EMed.gro 

cp $EQ1NVT .
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq1nvt.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq1nvt.mdp
sed -i 's/REFT/330 330 330/g' ./eq1nvt.mdp

gmx grompp -f eq1nvt.mdp -o SYSTEM-EQ1NVT.tpr -n SYSTEM.ndx -p SYSTEM-POSRES.top -maxwarn 10 -c SYSTEM-EMed.gro

$MDRUN  -deffnm SYSTEM-EQ1NVT -v -stepout 1000 -maxh 72 

cp $EQ1NPT .
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq1npt.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq1npt.mdp
sed -i 's/REFT/330 330 330/g' ./eq1npt.mdp

gmx grompp -f eq1npt.mdp -o SYSTEM-EQ1NPT.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ1NVT.cpt 

$MDRUN  -deffnm SYSTEM-EQ1NPT -v -stepout 1000 -maxh 72 

cp $EQ2 .
sed -i 's/GROUPS/DPPC SMA SOL/g' ./eq2.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./eq2.mdp
sed -i 's/REFT/330 330 330/g' ./eq2.mdp

gmx grompp -f eq2.mdp -o SYSTEM-EQ2.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ1NPT.cpt 

$MDRUN  -deffnm SYSTEM-EQ2 -v -stepout 1000 -maxh 72 

cp $MD .
sed -i 's/GROUPS/DPPC SMA SOL/g' ./md.mdp
sed -i 's/TAUT/1.0 1.0 1.0/g' ./md.mdp
sed -i 's/REFT/330 330 330/g' ./md.mdp

gmx grompp -f md.mdp -o SYSTEM-MD.tpr -n SYSTEM.ndx -p SYSTEM.top -maxwarn 10 -c SYSTEM.gro  -t SYSTEM-EQ2.cpt 

$MDRUN  -deffnm SYSTEM-MD -v -stepout 1000 -maxh 72 

rm \#* -f





























