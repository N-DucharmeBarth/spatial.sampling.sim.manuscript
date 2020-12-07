#!/bin/bash

pwd
ls -l
echo $PATH
set

# Upack everything from the tar file
tar -xzmf Start.tar.gz
tar -xzmf Background.tar.gz

# Untar the R executable
tar -xzf r361port.tar.gz

# dos2unix the script
dos2unix condor_run.Q.r
dos2unix rm_except
chmod 700 rm_except

# Export path variables
export PATH=$(pwd)/r361port/bin:$PATH
export R_ROOT_DIR=$(pwd)/r361port
export RHOME=${R_ROOT_DIR}

echo $PATH
echo $R_ROOT_DIR
echo $RHOME

# tests
# echo "Check"
# echo "What does 2*2 equal?"
# which R
# R RHOME
# Rscript -e 2*2

# Run the Rscript
Rscript --no-save --no-restore condor_run.Q.r > /dev/null 2>&1
# Rscript --no-save --no-restore condor_run.Q.r

./rm_except vast_list.RData
ls -l
