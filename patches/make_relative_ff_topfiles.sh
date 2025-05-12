#!/usr/bin/bash

for i in ./patch_*
do
	cd $i
	cp system.top system_rel.top
	sed -i 's|/mnt/3a29b482-dac1-4563-be89-d63ad92354e9/plasma_membrane/equil_chol/all_atom/d-erytho/backmapping-pipeline/patches|..|g' system_rel.top
	cd ..
done
