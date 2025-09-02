#!/bin/bash

path_output="./log/"


array=($(seq 0 1 9))

for index in "${array[@]}"; do

       echo $index
       sbatch --export=index=${index} ./GI_PPI_TotalOverlapSig_Run.sh

done
