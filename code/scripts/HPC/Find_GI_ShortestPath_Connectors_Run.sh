#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/home/xaviercg/scratch/log/%j.out
#SBATCH --export=index


# Ensure index is set before using it
if [ -z "$index" ]; then
    echo "Error: index is not set. Exiting."
    exit 1
fi

multiplier=$index

source ~/scratch/Yeast_Epistasis_Energetics_venv/bin/activate

echo 'TaskNumber'
echo $multiplier

# Generate an array from 0 to 7. The python script will be run 8 times in parallel.
array=($(seq 0 1 7))

# Generate the arguments for the Python script
args=()
for numerator in "${array[@]}"; do
    arg=$((multiplier * 8 + numerator))
    args+=("$arg")
done

# Function to run the Python script
run_python_script() {
    python ~/scratch/scripts/HPC_Find_GI_ShortestPath_Connectors.py "$1"
}

# Export the function to make it available to GNU Parallel
export -f run_python_script

# Run the Python script in parallel
parallel -j 8 run_python_script ::: "${args[@]}"

deactivate
