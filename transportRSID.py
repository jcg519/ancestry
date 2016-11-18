#/usr/bin/python
#!-*-coding:utf-8 -*-

import os,sys
import re
import math
import linecache 
from optparse import OptionParser


def parseCommand():
		usage = "usage: ./transportRSID.py -1 TM_merge_order_chr8.bim -o output"
		version = "%prog 1.0"
		parser = OptionParser(usage = usage, version = version)
		parser.add_option("-1", "--input1", dest = "input1", help = "the 23andme txt file")
		parser.add_option("-p", "--tped", dest = "tped", help = "the tped file")
		return parser.parse_args()

def getfile(input1,tped):
	count = len(open(input1).readlines())
	for i in range(1,count):
		fst_line=linecache.getline(input1,i)
		sec_line=linecache.getline(input1,i+1)
		fst=re.split('\t',fst_line)
		fst_rs=fst[1]
		sec=re.split('\t',sec_line)
		sec_rs=sec[1]
		fst_compare=fst[3]+fst[4]+fst[5]
		sec_compare=sec[3]+sec[4]+sec[5]
		if fst_compare==sec_compare:
			print fst_rs,sec_rs


if __name__=="__main___":
	time1=time.time()
	(options, args) = parseCommand()
	if options.input1 == None:
		sys.exit(-1)
	getfile(options.input1,options.output)
	time2=time.time()
	print 'Time used:' + str(time2-time1)

