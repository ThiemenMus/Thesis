qsub -pe serial $1 -l h_vmem=$2G -e logs/bc3.e -o logs/bc3.o bc3.txt