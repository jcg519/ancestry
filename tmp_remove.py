#/usr/bin/python
#!-*- coding:utf-8 -*-

import os,sys
import re
import linecache


count = len(open('TM_merge_order_chr9.bim').readlines())
f=open('TM_flib_chr9.tped','r').read()
for i in range(1,count):
	fst_line=linecache.getline('TM_merge_order_chr9.bim',i).strip()
	sec_line=linecache.getline('TM_merge_order_chr9.bim',i+1).strip()
	fst=re.split('\t',fst_line)
	fst_rs=fst[1]
	sec=re.split('\t',sec_line)
	sec_rs=sec[1]
	fst_compare=fst[3]+fst[4]+fst[5]
	sec_compare=sec[3]+sec[4]+sec[5]
	if fst_compare==sec_compare:
		print fst_rs,sec_rs
#		f=re.sub(sec_rs,fst_rs,f)
#		print f		

			
