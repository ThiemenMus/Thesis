qsub -pe serial $1 -l h_vmem=$2G -e logs/sergio.e -o logs/sergio.o Sergio.PBS
