#!/bin/csh
#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -r y

# The MATLABPATH variable is set in the Matlab script to add additional
# directories to the internal search paths.

#setenv MATLABPATH directory_path_to_your_files.m:other_user_contrib_directory_path
module load matlab

cd ..

matlab -singleCompThread -nodisplay -nosplash < testing.m
