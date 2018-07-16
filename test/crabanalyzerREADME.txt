To run the analyzer over multiple ntuples using a CRAB job:
1. Fill VLQAnalysis_Ntuples.txt with paths to ntuples, using AAA conventions.
2. Make certain you have the FrameworkJobReport.xml, myscript.sh, VLQAnalyzer.py, and crabConfig.py files all in the same directory.
3. Run: crab submit crabConfig.py

FrameworkJobReport.xml - Necessary file to run a CRAB job that doesn't run any cmsRun command.

myscript.sh - This shell script executes the proper commands to run the analyzer over the ntuples. May need to change the number of "if $1 -eq ?" lines to match the number of ntuples you are running over.

VLQAnalyzer.py - The analyzer code that we are executing.

crabConfig.py - The config file that creates and organizes this crab job. You should change the output directory to some storage site that you have writing access to, and you should change the totalunits to match the number of ntuples you are running over. The config file will pass each unit to the shell script as the first argument (in the my code, it passes numbers 1-29 to the shell script to separate the jobs).
