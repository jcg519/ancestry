import os,sys
import re

aa=open('some','r')
for line in aa.readlines():
	tmp=re.split("\s+",line.strip())
	rsid=tmp[1]
	alt=tmp[4]
	ref=tmp[5]
	gt=tmp[7:]
	hom=ref+ref
	het1=alt+ref
	het2=ref+alt
	ho =alt+alt
	ch ='00'
	rep = ['0 0 1' if x ==hom or x==ch else x for x in gt]
	total=['0 1 0' if x ==het1 else x for x in rep]
	hh=['0 1 0' if x ==het2 else x for x in total]
	end =['1 0 0' if x ==ho else x for x in hh]
	strlist=' '.join(end)
	print tmp[0],tmp[1],tmp[3],tmp[4],tmp[5],strlist

