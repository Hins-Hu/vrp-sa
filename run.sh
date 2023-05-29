#/bin/bash

# Define the result folder
output_folder="Results-P/"

# Define the set of inputs
inputs=("P-n101-k4")


# Iterate over the inputs
for input in "${inputs[@]}"
do

    path="${output_folder}${input}"

    if [ -d "$path" ]; then
        rm $path/*
    else
        mkdir -p "$path"
    fi

    output_file="${output_folder}${input}/logging.log"
    nohup python main.py "$input" "$path" > "$output_file" &
done


