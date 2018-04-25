#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
import sys
#print "This is the name of the script: ", sys.argv[0]
#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

#ccr = ' resubmit '
#ccs = ' status '
#cck = ' kill '

if( len(sys.argv) != 3 ): 

	print "Wrong number of arguments"
	sys.exit()

cc = sys.argv[2]

#thedir = "80X_trees_Sig_cleaned_08Aug17/"
#thedir = "80X_trees_Test2_25Aug17/"
thedir = sys.argv[1]

os.system( "ls "+thedir+" > filelist.txt" )
infile = open( "filelist.txt", "r" )
for line in infile:
	print ">>>>>>>>>>>>>>>>>> crab" + cc + " : "+thedir+line
	os.system( "crab "  + cc + " " + thedir + line )	
