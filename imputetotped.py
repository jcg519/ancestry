#/usr/bin/python
#!-*- coding:utf-8 -*-

#--author--=='cjiang'

##############trans the impute  file to tped file ###### 

############## the impute file ##########
'''
--- rs149630736:45000028:T:C 45000028 T C 1 0 0 1 0 0 1 0 0 1 0 0
--- rs146905936:45000294:G:C 45000294 G C 1 0 0 1 0 0 1 0 0 1 0 0 
10 rs1253750 45000220 0 G 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 
 '''
#########################################

import os,sys
import re
import time
import optparse



def DNA_complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()
if __name__=='__main__':
	usage = "imputetotped [ -i <genfule> -o <outdir>]"
	opter = optparse.OptionParser(usage)
	opter.add_option("-i", "--file", dest="input", help="all gen file")
	opter.add_option("-o", "--output", dest="output", help="to tped file")
	opter.add_option("-r", "--ref", dest="ref", default="reference.bim"  ,help="the reference file")
	opt, args = opter.parse_args()
	input = opt.input
	output = opt.output
	ref = opt.ref
	a=open(input,'r')
	b=open(ref,'r')
	c=open(output,'w')
	chr={}
	chgeno={}
	for line in b.read():
		tmp=re.split("\t",line)
		chr[tmp[1]]=tmp[0]
		chgeno[tmp[1]]=[tmp[4],tmp[5].strip()]

	count=0
	old='0'
	with open(input) as f:
		for lines in f.readlines():
			res=re.split("\s+",lines)
			rsid=re.split(":",res[1])[0]
			new=rsid
			pos=res[2]
			ref=res[3]
			alt=res[4]
			genty=res[5:]
			chrom = chr.get(rsid)
			if chrom:
				new=rsid
				if new !=old:
					lh=len(genty)/3
					genetype=[]
					for i in range(0,lh):
						i1=3*i
						i2=3*i+1
						i3=3*i+2
						reflist = chgeno.get(rsid)
						if ref in reflist or alt in reflist:
							if genty[i1]>=0.9 or genty[i2]>=0.9 or genty[i3]>=0.9:
								if genty[i1]>=genty[i2] and genty[i1]>=genty[i3]:
									gene=ref+" "+ref
									genetype.append(gene)
								if genty[i2] > genty[i1] and genty[i2] >= genty[i3]:
									gene=ref+" "+alt
									genetype.append(gene)
								if genty[i3] > genty[i1] and genty[i3] >genty[i2]:
									gene=alt+" "+alt
									genetype.append(gene)
							else:
								pass
						elif DNA_complement(ref) in reflist and DNA_complement(alt) in reflist:
							if genty[i1]>=0.9 or genty[i2]>=0.9 or genty[i3]>=0.9:
						
								if genty[i1]>=genty[i2] and genty[i1]>=genty[i3]:
									gene=DNA_complement(ref)+" "+DNA_complement(ref)
									genetype.append(gene)
								if genty[i2] > genty[i1] and genty[i2] >= genty[i3]:
									gene=DNA_complement(ref)+" "+DNA_complement(alt)
									genetype.append(gene)
								if genty[i3] > genty[i1] and genty[i3] >genty[i2]:
									gene=DNA_complement(alt)+" "+DNA_complement(alt)
									genetype.append(gene)
							else:
								pass

					gtgen=" ".join(genetype)
					if len(gtgen)>1:
						lines=chrom+" "+rsid+" "+'0'+" "+pos+" "+gtgen+"\n"
						c.writelines(lines)
	old=new
	count+=1
	
	a.close()
	b.close()
