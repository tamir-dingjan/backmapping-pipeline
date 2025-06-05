#!/bin/bash

touch eq.logging
echo "Starting patch equilibration" > eq.logging

for patchdir in /home/labs/futerlab/tamird/equil_chol_all_atom/3keto/patches/patch_*/
do
	cd $patchdir
	if test -f "md0.tpr"; then
		echo "$patchdir is equilibrated" >> ../eq.logging
	else
		cat > eq.submit <<EOF
#BSUB -L /bin/bash
#BSUB -q long-gpu
#BSUB -J eq
#BSUB -oo eq_output
#BSUB -eo eq_error
#BSUB -gpu "num=1:j_exclusive=yes"
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "affinity[thread(20, same=socket)]"
#BSUB -W 5750
#BSUB -cwd $patchdir
EOF
		cat >> eq.submit <<'EOF'
module spider GROMACS > gmxmodules |& tee
module load $(grep -m 1 CUDA gmxmodules | tail -n1)

if [ ! -f nvt.gro ]; then
  gmx grompp -f ../nvt.mdp -c min.gro -p system_rel.top -n index.ndx -o nvt.tpr
  gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 20
fi

if [ ! -f npt.gro ]; then
  gmx grompp -f ../npt.mdp -c nvt.gro -p system_rel.top -n index.ndx -o npt.tpr
  gmx mdrun -deffnm npt -ntmpi 1 -ntomp 20
fi

# When automatically selecting a Gromacs module, cannot rely on the TPR here 
# being readable by the production job
#if [ ! -f md0.tpr ]; then
#  gmx grompp -f ../md0.mdp -c npt.gro -p system_rel.top -n index.ndx -o md0.tpr
#fi

EOF

	fi

	bsub < eq.submit

	cd ..
done

