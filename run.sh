#/bin/bash

# Define the set of inputs
# inputs=("P-n20-k2" "P-n21-k2" "P-n22-k2" "P-n22-k8" "P-n23-k8" "P-n40-k5" "P-n45-k5" "P-n50-k7" "P-n50-k8" 
#         "P-n50-k10")

inputs=("X-n101-k25")
# Iterate over the inputs
for input in "${inputs[@]}"
do
    output_file="${input}.txt"
    nohup python main.py "$input" > "$output_file" &
done
