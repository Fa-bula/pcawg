#!/bin/bash
# Script to perform all needed calculations and wait for it finish

ROOT=../			# Dir with raw_mutations.txt and enrichment.txt
LOGS=$ROOT/logs/		# qsub .log files directory

rm $LOGS/*			# Remove old .log files
clear				# Clear terminal screen
date > $LOGS/exec_time	# Write start time at file

qsub -t 23-24 -tc 4 -sync y -cwd -S /usr/bin/python -e $LOGS -o $LOGS main.py

echo "Job has been finished!"
date >> $LOGS/exec_time		# Append end time to file
cat $LOGS/exec_time		# Show running start and end times
