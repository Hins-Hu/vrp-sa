#/bin/bash

# Define the set of instances
instances=("P-n16-k8" "P-n19-k2" "P-n20-k2" "P-n21-k2" "P-n22-k2" "P-n22-k8" "P-n23-k8"
        "P-n40-k5" "P-n45-k5" "P-n50-k7" "P-n50-k8" "P-n50-k10" "P-n51-k10" 
        "P-n55-k7" "P-n55-k10" "P-n55-k15" "P-n60-k10" "P-n60-k15" "P-n65-k10" 
        "P-n70-k10" "P-n76-k4" "P-n76-k5" "P-n101-k4") 

# instances=("P-n51-k10")

# Define the set of t_max_factor
# t_max_factors=("1.2") 

# Define the set of budget_factor
# budget_factors=("1")

# Define the set of grid_size
# grid_size=("2" "3" "7" "8")

# Cost adjustment
# eta_1=("0.3" "0.35" "0.4" "0.45" "0.5" "0.55" "0.6" "0.65" "0.7" "0.75" "0.8" "0.85" "0.9")
# eta_2=("1.7" "1.65" "1.6" "1.55" "1.5" "1.45" "1.4" "1.35" "1.3" "1.25" "1.2" "1.15" "1.1")
eta_1=("0.5" "0.6" "0.7" "0.8")
eta_2=("1.2" "1.3" "1.4" "1.5")


for e_1 in "${eta_1[@]}"
do
    for e_2 in "${eta_2[@]}"
    do 

# for g in "${grid_size[@]}"
# do
        # Define the result folder
        output_folder="Results-cost/Results-${e_1}-${e_2}"

        # Iterate over the inputs
        for i in "${instances[@]}"
        do

            path="${output_folder}/${i}"

            if [ -d "$path" ]; then
                rm $path/*
            else
                mkdir -p "$path"
            fi

            output_file="${output_folder}/${i}/logging.log"
            nohup python main.py "$i" "-p" "$path" "-e1" "$e_1" "-e2" "$e_2" "-g" "5"> "$output_file" &
        done
# done
    done
done