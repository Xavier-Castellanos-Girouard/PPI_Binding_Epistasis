#!/bin/bash

path_output="./log/"


array=($(seq 0 1 124))

for index in "${array[@]}"; do

       echo $index
       sbatch --export=index=${index} ./Find_GI_ShortestPath_Connectors_Run.sh

done
