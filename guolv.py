import os,sys
import re
import time


with open('test.guolv') as f:
    for line in f:
        tmp=re.split("\s+",line.strip())
        gt=tmp[5:]
        ls=len(gt)/3
        mostmax=[]
        for i in range(0,ls):
            ll=max(gt[3*i],gt[3*i+1],gt[3*i+2])
            mostmax.append(ll)
#        if len(mostmax)==95:
        if float(min(mostmax))>0.9:
            print line
        else:
            continue








