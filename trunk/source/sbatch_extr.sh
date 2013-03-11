#!/bin/bash

# This is the template sbatch script for Lab2 at BB2490 spring 2013

# Run the script like this:
# sbatch <OPTIONS> sbatch_template.sh
# Where the OPTIONS are all of:
# –A g2012230 –t 2:00:00 –p node –n 8 --res=g2012230_1 –o my_job_name.out
# (you may change the job name "my_job_name.out" to anything you like).

# Any line that starts with '# ' (dash followed by space) is a comment and will not be executed

# Don't forget to add the appropriate modules! See THE MODULE SYSTEM section in the instructions.

# cd to the directory where the data is (of course, use your own user name):
cd /home/nguyenh/rt-errors-prediction/source/

# Then list the commands that you would like to run. 
# If you'd like to add more jobs, enter them on the following lines. 
# To help you out, the first command of today's lab is provided:
time R --slave -f extr.r >> ../results/2mer/log 2>&1
