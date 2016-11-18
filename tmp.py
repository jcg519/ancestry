import os,sys
import re

aa=open('all.tped','r')
for line in aa.readlines():
	res=re.split("\s+",line)
	gt=res[4:]
	subject="".join(gt)
	result = re.sub(r"(?<=\w)(?=(?:\w\w)+$)", " ", subject)
	print result



