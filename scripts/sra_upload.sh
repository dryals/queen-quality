#!/bin/bash

# FILENAME: sra_upload.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=cpu
#SBATCH --job-name sra
#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/sra_upload.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/sra_upload.out

echo "moving..."
cd /scratch/negishi/dryals/queen-quality/fastq/SRA-2022-2026

echo "lftp..."

lftp -u subftp,Eg8FrojfichDymhadfok ftp-private.ncbi.nlm.nih.gov << EOF

cd uploads/dylan.k.ryals_gmail.com_xQqWLtyf/SUB123456_related_data_2
mput *

bye
EOF


echo "DONE"


