#!/bin/bash

path_output="./log/"


array=($(seq 0 1 9))

for index in "${array[@]}"; do

       echo $index
       sbatch --export=index=${index} ./Identify_Connector_Overlap_Run.sh

done
