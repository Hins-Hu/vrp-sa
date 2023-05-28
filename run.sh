#/bin/bash

# Define the set of inputs
inputs=("A-n33-k6")

# Iterate over the inputs
for input in "${inputs[@]}"
do
    output_file="${input}.txt"
    nohup python main.py "$input" > "$output_file" &
done


