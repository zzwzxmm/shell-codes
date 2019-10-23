#!/bin/bash
source /home/xuemm/software/gromacs-5.0.7/bin/GMXRC
GMXGMP='gmx grompp'
MDRUN='gmx mdrun'
array_ele=(`seq -0.010 0.001 0.010`)
CP=/bin/cp
RM=/bin/rm
MV=/bin/mv

for ele in ${array_ele[@]}
do
DIR="E_$ele"
if [ -e $DIR ]; then
cd $DIR
else
mkdir $DIR && cd $DIR
fi
ELE=$ele
$CP ../md.mdp md.mdp
$CP ../gmx.qsub gmx.qsub
$CP -R ../charmm.ff/ .
sed -i "s/0 0.5 0/1 $ele 0/g" md.mdp
#n_p=`cat /proc/cpuinfo | grep processor | wc -l`
#n_omp=$((n_p/2))
sed -i "s/JOBNAME/E_$ele/g" gmx.qsub
#sed -i "s/NUMBEROFPROCESSORS/$n_p/g" gmx.qsub
#sed -i "s/NUMBEROFOMPTHREADS/$n_omp/g" gmx.qsub
$GMXGMP -f md.mdp -c ../em.gro -n ../system.ndx -p ../system.top -maxwarn 10 -o md.tpr
qsub gmx.qsub
cd ..
done


#for ele in ${array_ele[@]}
#do
#DIR="E_$ele"
#cd $DIR
#echo "mdrun -deffnm md -v &> md.log"
#$MDRUN -deffnm md -v -maxh 240 &> md.log
#cd ..
#done

rm -f \#mdout* 
