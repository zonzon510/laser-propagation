#!/bin/bash

# calculate the driving laser field
laser_propagation/build/main.out ./freespace.input

# run hhgmax with the output
mkdir -p driving_field
matlab -nodesktop -r calculate_harmonics

