#!/bin/bash

#a bash script to generate the flux of certain molecule groups
#based on gromacs-4.5.6 and github tool g_count/g_flux/g_zcoord

source /home/xuemm/software/gromacs-4.5.6/bin/GMXRC
for i in `seq -0.5 0.1 0.5`
do
DIR=E_$i
cd $DIR
grompp -f md.mdp -p ../top.top -n ../system.ndx -c ../eq.gro -maxwarn 10 -o md_g456.tpr
echo -en "5\n" | g_flux -f md.xtc -s md_g456.tpr -n ../system.ndx -o flux_cl.xvg -R 5 -z1 3 -z2 8 -b 10000 -e 50000
echo -en "4\n" | g_flux -f md.xtc -s md_g456.tpr -n ../system.ndx -o flux_k.xvg -R 5 -z1 3 -z2 8 -b 10000 -e 50000
tail -n 1 flux_cl.xvg >> ../flux_cl_0.1.dat
tail -n 1 flux_k.xvg >> ../flux_k_0.1.dat
cd ..
done


for i in `seq -0.5 0.1 0.5`
do
read -a line
echo $i ${line[6]}
done < flux_cl_0.1.dat > current_cl_0.1.xvg

for i in `seq -0.5 0.1 0.5`
do
read -a line
echo $i ${line[6]}
done < flux_k_0.1.dat > current_k_0.1.xvg
