#!/bin/bash

touch md0.logging
echo "Starting patch simulation md0" > md0.logging

for patchdir in /home/labs/futerlab/tamird/equil_chol_all_atom/d-erythro/patches/patch_*/
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
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "affinity[thread(10, same=socket)]"
#BSUB -W 5750
#BSUB -cwd $patchdir

module load GROMACS/2023.1-foss-2022a-CUDA-11.7.0

if [ ! -f md0.cpt ]; then
  gmx mdrun -deffnm md0 -ntmpi 2 -ntomp 10 -cpo md0.cpt  
else
  gmx mdrun -deffnm md0 -ntmpi 2 -ntomp 10 -cpi md0.cpt -cpo md0.cpt
fi

EOF

	fi

	bsub < md0.submit

	cd ..
done

