qsub -pe serial $1 -l h_vmem=$2G -e logs/genie.e -o logs/genie.o genie.txt