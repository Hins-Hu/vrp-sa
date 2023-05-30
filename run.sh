#/bin/bash

# Define the set of instances
instances=("P-n16-k8" "P-n19-k2" "P-n20-k2" "P-n21-k2" "P-n22-k2" "P-n23-k8"
        "P-n40-k5" "P-n45-k5" "P-n50-k7" "P-n50-k8" "P-n50-k10" "P-n51-k10" 
        "P-n55-k7" "P-n55-k10" "P-n55-k15" "P-n60-k10" "P-n60-k15" "P-n65-k10" 
        "P-n70-k10" "P-n76-k4" "P-n76-k5" "P-n101-k4") 
# instances=("P-n50-k7")

# Define the set of t_max_factor
t_max_factors=("1" "1.2" "1.5")

# Define the set of budget_factor
budget_factors=("1" "2" "3")



for t_f in "${t_max_factors[@]}"
do 
    for b_f in "${budget_factors[@]}"
    do


        # Define the result folder
        output_folder="Results-${t_f}-${b_f}"


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
            nohup python main.py "$i" "$path" "-t" "$t_f" "-bf" "$b_f"> "$output_file" &
        done


    done
done


