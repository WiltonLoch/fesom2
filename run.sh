#! /bin/sh
#

user=k202167

# Remove restart files to avoid restart run
rm results/*restart*

# Set the fesom clock to the initialization value
cp fesom.clock results/fesom.clock

# Enter the work folder to avoid creation of files in the root folder
cd work

export NVCOMPILER_ACC_NOTIFY=3

if [ "$1" = "profile" ]
then
    export NSYS_PROFILE=1
fi

sbatch job_levante
squeue -i 1 -u $user
