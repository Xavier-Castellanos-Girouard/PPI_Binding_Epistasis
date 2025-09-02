#!/bin/bash

#Clustal v2.1

# Execute the script in results/IntermediateFiles directory

for f in Paralog_proteinSeq_ClustalW_input/*.fasta ; do
	echo "$f"
	
	prefix="Paralog_proteinSeq_ClustalW_input/"
	suffix=".fasta"
	output_name=${f#"$prefix"}
	output_name=${output_name%"$suffix"}
	output_dir="ClustalW_Output/${output_name}.aln"
	
	#echo "${output_dir}"
	clustalw -INFILE=$f -ALIGN -outfile=$output_dir -type=PROTEIN -seqnos=ON -pwmatrix=PAM -quiet
done
