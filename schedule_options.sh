#!/bash/bin

# Make directory for output files
mkdir -p output_options

# Run program
python3 -W"ignore" scheduling_options.py

# Print schedule
echo "Schedule:"
cat sched_realta_out.txt

# Remove the temporary files
rm sched_iLiSA_temp.txt
rm sched_realta_temp.txt

# Move the created files to the output folder
mv sched_iLiSA_out.txt ./output_options
mv sched_realta_out.txt ./output_options
mv optimum.png ./output_options
mv viewing.png ./output_options
