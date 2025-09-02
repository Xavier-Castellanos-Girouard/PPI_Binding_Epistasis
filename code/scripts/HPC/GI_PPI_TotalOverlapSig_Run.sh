#!/bin/bash

#SBATCH --time=2:30:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/xaviercg/scratch/log/%j.out
#SBATCH --export=index

# Ensure index is set before using it
if [ -z "$index" ]; then
    echo "Error: index is not set. Exiting."
    exit 1
fi

multiplier=$index

source ~/scratch/Yeast_Epistasis_Energetics_venv/bin/activate

echo 'Task Number'
echo $multiplier

# Generate an array from 0 to 999. The python script will be run 10000 times.
array=($(seq 0 1 999))

# Generate the arguments for the Python script
args=()
for numerator in "${array[@]}"; do
    arg=$((multiplier * 1000 + numerator))
    args+=("$arg")
done

# Function to run the Python script
run_python_script() {
    python ~/scratch/scripts/HPC_Module_Overlap_Analysis.py "$1"
}

# Export the function to run the Python script
export -f run_python_script

# Run the Python script in parallel
parallel -j 8 run_python_script ::: "${args[@]}"

deactivate
