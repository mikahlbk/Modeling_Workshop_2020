#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-2:00:00
#SBATCH --output=my1c_40.stdout
#SBATCH --job-name="1c_40"
#SBATCH -p short 
export OMP_NUM_THREADS=12
mkdir /bigdata/wchenlab/shared/Workshop/Animate_Cyt_1c_40
mkdir /bigdata/wchenlab/shared/Workshop/Nematic_test_1c_40
mkdir /bigdata/wchenlab/shared/Workshop/Locations_test_1c_40
mkdir /bigdata/wchenlab/shared/Workshop/Animate_No_Cyt_1c_40
mkdir /bigdata/wchenlab/shared/Workshop/Cell_Data_1c_40
mkdir /bigdata/wchenlab/shared/Workshop/Tissue_Data_1c_40
./program /bigdata/wchenlab/shared/Workshop/Animate_Cyt_1c_40 /bigdata/wchenlab/shared/Workshop/Locations_test_1c_40 /bigdata/wchenlab/shared/Workshop/Nematic_test_1c_40 /bigdata/wchenlab/shared/Workshop/Animate_No_Cyt_1c_40 /bigdata/wchenlab/shared/Workshop/Cell_Data_1c_40 /bigdata/wchenlab/shared/Workshop/Tissue_Data_1c_40 -div 3 -Bending 25
