#! /bin/bash
#
# Run GETM on multiple CPUs
#
# To run this script in the background with nohup, use:
# $ nohup ./run.sh &
#
# by Markus Reinert (IOW), April 2022 to March 2023

# Set the path to the GETM code directory
GETMDIR=$HOME/tools/getm/code/ 

# Choose the number of CPU cores to use
nCPU=4

echo "Creating input files"
python3 Make_fjord_setup_322.py

echo "Creating namelist files"
editscenario --schemadir $GETMDIR/schemas -e nml . fjord_322.xml

echo "Starting GETM ..."
time mpirun -n $nCPU bin/getm

echo "Merging snapshots ..."
ncmerge store/out_snapshot.00*.nc store/out_snapshot.nc && rm store/out_snapshot.00*.nc
echo "Merging means ..."
ncmerge store/out_mean.00*.nc store/out_mean.nc && rm store/out_mean.00*.nc

echo "Finished!"
