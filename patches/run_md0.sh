#!/bin/bash

touch md0.logging
echo "Starting patch simulation md0" > md0.logging

for patchdir in /home/labs/futerlab/tamird/equil_chol_all_atom/3keto/patches/patch_*/
do
	cd $patchdir
	if test -f "md0.gro"; then
		echo "$patchdir is complete" >> ../md0.logging
	else
		cat > md0.submit <<EOF
#BSUB -L /bin/bash
#BSUB -q long-gpu
#BSUB -J md0
#BSUB -oo md0_output
#BSUB -eo md0_error
#BSUB -gpu "num=1:j_exclusive=yes"
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "affinity[thread(20, same=socket)]"
#BSUB -W 5750
#BSUB -cwd $patchdir

module spider GROMACS >> gmxmodules |& tee
module load $(grep -m 1 CUDA gmxmodules | tail -n1)

if [ ! -f md0.cpt ]; then
  gmx grompp -f ../md0.mdp -c npt.gro -p system_rel.top -n index.ndx -o md0.tpr
  gmx mdrun -deffnm md0 -ntmpi 1 -ntomp 20 -cpo md0.cpt  
else
  gmx mdrun -deffnm md0 -ntmpi 1 -ntomp 20 -cpi md0.cpt -cpo md0.cpt
fi

EOF

	fi

	bsub < md0.submit

	cd ..
done

