#!/bin/bash

# Make directory for output files
mkdir -p output

# Make the schedule
python3 -W"ignore" sched.py $1 $2 $3 $4

# Remove the quotation marks
sed -i 's/"//g' sched_realta.txt
sed -i 's/"//g' sched_iLiSA.txt

# Print schedule
echo "Schedule:"
cat sched_realta.txt

# Move the created files to the output folder
mv sched_iLiSA.txt ./output
mv sched_realta.txt ./output
