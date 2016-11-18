#!/usr/bin/env python
#coding:utf-8
__author__ = 'similarface'
import sys
import os
import re
import math
import optparse
from collections import defaultdict

class SplitGenFile:
    def splitChromGen(self,genfile):
        '''
        12      rs13273 133498344       G       G       1       0       0
        14      rs7157079       95327349        A       A       1       0       0
        :param genfile:
        :return:
        '''
        chrdict=defaultdict(list)
        for line in open(genfile,'r'):
            lines=re.split('\s+',line.strip())
            chrdict[lines[0]].append(lines)
        return chrdict

    def split5MChromGen(self,chrdict):
        result=defaultdict(list)
        chrlen={'1':249250621,'2':243199373,'3':198022430,'4':191154276,'5':180915260,'6':171115067,'7':159138663,'8':146364022,'9':141213431,'10':135534747,'11':135006516,'12':133851895,'13':115169878,'14':107349540,'15':102531392,'16':90354753,'17':81195210,'18':78077248,'19':59128983,'20':63025520,'21':48129895,'22':51304566}
        for i in range(1,23):
            data=chrdict[str(i)]
            widthlen=math.ceil(chrlen[str(i)]/5000000.0)
            for j in range(1,int(math.ceil(chrlen[str(i)]/5000000.0))+1):
                for item in data:
                    if (j-1)*5000000<=int(item[2])<=j*5000000:
                        result["chr"+str(i)+"_"+str(j)].append(item)
        return result

    def getSplitFile(self,filesavedict,genfile):
        if not os.path.exists(filesavedict):
            os.system("mkdir -p "+filesavedict)
        result=self.split5MChromGen(self.splitChromGen(genfile))
        for k,v in result.items():
            filename=os.path.join(filesavedict,k+".gen")
            fileoper=open(filename,'w')
            try:
                for iten in v:
                    fileoper.writelines('\t'.join(iten)+'\n')
                fileoper.close()
            except Exception,e:
                fileoper.close()

if __name__ == "__main__":
    usage = "splitgen [ -f <genfule> -d <outdir>]"
    opter = optparse.OptionParser(usage)
    opter.add_option("-f", "--file", action="store", dest="genfile", help=u"生成的gen文件")
    opter.add_option("-d", "--dir", action="store", dest="outdir", help=u"切割输出目录")
    opt, args = opter.parse_args()
    genfile = opt.genfile
    outdir = opt.outdir
    if genfile != None and outdir != "" and outdir!=None and outdir!="":
         spgen=SplitGenFile()
         spgen.getSplitFile(outdir,genfile)
    else:
        opter.print_help()
        sys.exit(0)
