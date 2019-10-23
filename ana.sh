#!/bin/bash
#source gmx
source /home/xuemm/software/gromacs-5.0.7/bin/GMXRC
#analyze the current using gmx current
for i in `seq -0.010 0.001 0.010`
do 
cd E_$i 
echo -en "9\nq\n" | gmx current -s md.tpr -n ../system.ndx -f md.trr -nojump -o g_cur.xvg 
cd ..
done

for i in `seq -0.010 0.001 0.010`
do 
cd E_$i 
gmx analyze -f mj.xvg &>> ../mj-ana.log
cd .. 
done

cat mj-ana.log | grep SS3 > mj-ana.xvg

for i in `seq -0.010 0.001 0.010`
do
read -a line
echo $i ${line[1]} ${line[3]}
done < mj-ana.xvg > mj-ana-cut.xvg
