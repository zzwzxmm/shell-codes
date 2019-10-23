#!/bin/bash
#source gmx
#source /home/xuemm/local/gromacs/bin/GMXRC
#analyze the current using gmx current
for i in `seq -0.10 0.01 0.10`
do 
cd E_$i 
#echo -en "7\nq\n" | gmx current -s md.tpr -n ../system.ndx -f md.trr -nojump -o g_cur.xvg 
runvmd ../current.tcl
cd ..
done

for i in `seq -0.10 0.01 0.10`
do 
cd E_$i 
#gmx analyze -f mj.xvg &>> ../mj-ana.log
tail -n 1 >> ../current.xvg
cd .. 
done

