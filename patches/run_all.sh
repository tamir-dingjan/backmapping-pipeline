#!/bin/bash

touch sim.logging
echo "Starting patch simulation equilibration and production" > sim.logging

for patchdir in /home/labs/futerlab/tamird/equil_chol_40molpc_all_atom_correct_sn1_sn2/d-erythro/patches/patch_*/
do
    cd $patchdir
    if test -f "md0.gro"; then
        echo "$patchdir has finished production simulation" >> ../sim.logging
    else
        cat > run.submit <<EOF
#BSUB -L /bin/bash
#BSUB -q long-gpu
#BSUB -J sim
#BSUB -oo sim_output
#BSUB -eo sim_error
#BSUB -gpu "num=1:j_exclusive=yes"
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "affinity[thread(20, same=socket)]"
#BSUB -W 5750
#BSUB -cwd $patchdir
EOF
        cat >> run.submit <<'EOF'
module spider GROMACS >> gmxmodules |& tee
module load $(grep -m 1 CUDA gmxmodules | tail -n1)

if [ ! -f nvt.gro ]; then
  gmx grompp -f ../nvt.mdp -c min.gro -p system_rel.top -n index.ndx -o nvt.tpr
  gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 20
fi

if [ ! -f npt.gro ]; then
  gmx grompp -f ../npt.mdp -c nvt.gro -p system_rel.top -n index.ndx -o npt.tpr
  gmx mdrun -deffnm npt -ntmpi 1 -ntomp 20
fi

if [ ! -f md0.gro ]; then
  gmx grompp -f ../md0.mdp -c npt.gro -p system_rel.top -n index.ndx -o md0.tpr
  gmx mdrun -deffnm md0 -ntmpi 1 -ntomp 20 -cpo md0.cpt 
fi

EOF

    fi

    bsub < run.submit

    cd ..
done
