#!/bin/bash

for patch in ./patch_*
do
	if [ -f $patch/md0.gro ]; then
		if grep -q False $patch/stereo.check; then
			echo $patch
		fi
	fi
done

			
