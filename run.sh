#!/bin/bash

# calculate the driving laser field
laser_propagation/build/main.out ./freespace.input

# run hhgmax with the output
matlab -nodesktop -r calculate_harmonics

