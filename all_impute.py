import os,sys
import splitgenfile
import re
import time
import all
from optparse import OptionParser
from multiprocessing import Process, Queue
import multiprocessing


def parseCommand():
		usage = "usage: ./impute.py -1 inputfile  -n name -o outputdir"
		version = "%prog 1.0"
		parser = OptionParser(usage = usage, version = version)
		parser.add_option("-1", "--inputfile", dest = "inputfile", help = "the split gen file dir")
		parser.add_option("-n", "--name", dest="name", help="the name file ")
		parser.add_option("-o", "--outputdir", dest = "outputdir", help = "the output file dir")
		return parser.parse_args()


def getgen(inputfile,name, output):
	os.system("/usr/bin/python /data4/AncestryAnalysis/tools/gen.py -i "+ inputfile + " -o " + name)
	os.system("/usr/bin/python /data4/AncestryAnalysis/tools/splitgenfile.py -f "+ name+ " -d "+ output)
	os.system("/usr/bin/python /data4/AncestryAnalysis/tools/all.py -o "+ output)
	return 1

if __name__=="__main__":
	time1 = time.time()
	(options, args) = parseCommand()
	getgen(options.inputfile, options.name, options.output)
	time2 = time.time()
	print "time_used:"+str(time2-time1)






