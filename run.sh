#! /bin/sh
#

user=k202167

# Remove restart files to avoid restart run
rm results/*restart*

# Set the fesom clock to the initialization value
cp fesom.clock results/fesom.clock

# Enter the work folder to avoid creation of files in the root folder
cd work
sbatch job_levante
timeout 20s squeue -i 1 -u $user
