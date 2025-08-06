#!/usr/bin/bash

for i in ./patch_*
do
	cd $i
	cp system.top system_rel.top
	sed -i 's|/media/tamir/Elements/simulation_data_backups/plasma_membrane/equil_chol_40molpc/all_atom_correct_sn1_sn2/d-erythro/backmapping-pipeline/patches|..|g' system_rel.top
	cd ..
done
