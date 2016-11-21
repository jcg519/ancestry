#!/usr/bin/python

import os,sys,re
import logging
import commands
from optparse import OptionParser

def parseCommand():
	usage = "usage: %prog <-1 input1> <-2 input2> <-o output>"
	version = "%prog 1.0"
	parser = OptionParser(usage = usage, version = version)
	parser.add_option("-1", "--input1", dest = "input1",
		help = "the first annovar vcf file")
	parser.add_option("-v", "--varq", dest = "varQuality",
		help = "set varQuality threshold")
	parser.add_option("-g", "--genome1000", dest = "genome1000",
		help = "set 1000 genome filter")
	parser.add_option("-o", "--output", dest = "output", default = "common.anno.csv",
		help = "the common annovar output csv file")
	return parser.parse_args()

if __name__ == "__main__":
	(options, args) = parseCommand()
	if options.input1 == None:
		print "input1 is None"
		sys.exit()
	if options.output == None:
		print "output is None"
		sys.exit()
	genome1000 = 0
	if options.genome1000:
		genome1000 = float(options.genome1000)
#		print genome1000
#		sys.exit()
	fd_w = open(options.output, "w")
	with open(options.input1, "r") as fd:
		for line in fd.readlines():
			if re.match("^#", line):
				fd_w.write(line)
			elif re.match(".*;snp138=\.", line):
				if options.genome1000:
					m = re.match(".*;1000g2014oct_all=(.*?);", line)
					if m:
						if m.group(1) == ".":
							fd_w.write(line)
						elif float(m.group(1)) <= genome1000:
							fd_w.write(line)
#			break
#	print "hello"
	fd.close()  
	fd_w.close()

