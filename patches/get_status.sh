#!/usr/bin/bash
watch -n 1 "tail -n 15 ./patch_*/md0.log | grep -e == -e Time -A1 -e Fatal -e Finished | grep -e == -e 0000 -e Fatal -e Finished"

