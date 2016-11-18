#coding:utf-8
__author__ = 'similarface'

import time
import math
import re
import os
import vcf
import sys
from collections import defaultdict
import MySQLdb

class VcfComvertToOtherFile:

    def __init__(self,vcffilename):
        self.vcffilename=vcffilename
        try:
            self.__db=MySQLdb.connect(host='192.168.30.252',port=3306,db='gendb',user='dna',passwd='dna',charset='utf8')
            self.__cursor = self.__db.cursor()
        except MySQLdb.Error, e:
          self.error_code = e.args[0]
          error_msg = 'MySQL error! ', e.args[0], e.args[1]
          print error_msg
    #
    # def __enter__(self):
    #     '''
    #     获取数据库连接 __enter__函数配合with 使用 自动调用 enter 和 exit
    #     :return:
    #     '''
    #     self.__db=MySQLdb.connect(host='192.168.30.252',port=3306,db='gendb',user='dna',passwd='dna',charset='utf8')
    #     self.__cursor = self.__db.cursor()
    #     return self
    #
    # def __exit__(self, type, value, traceback):
    #     self.__db.close()


    def vcfTo23andmeFile(self,outFile):
        pass

    def vcfToWegeneFile(self,outFile):
        pass

    def vcfToPMRAfile(self,outFile):
        pass

    def get23andmeFiledict(self,andmefile):
        '''
        :param andmefile: 输入23andme文件
        :return:key: chr:pos value:rsid 的字典
        '''
        regx=re.compile('\s+')
        andmefiledict={}
        with open(andmefile,'r') as oper:
            for line in oper:
                if line.startswith("#") or line.strip()=="":
                    pass
                else:
                    lines=re.split(regx,line.strip())
                    andmefiledict[lines[1]+":"+lines[2]]=lines[0]
        return andmefiledict

    def getwegeneFiledict(self,wegenefile):
        '''
        :param andmefile: 输入23andme文件
        :return:key: chr:pos value:rsid 的字典
        '''
        regx=re.compile('\s+')
        wegenefiledict={}
        with open(wegenefile,'r') as oper:
            for line in oper:
                if line.startswith("#") or line.strip()=="":
                    pass
                else:
                    lines=re.split(regx,line.strip())
                    wegenefiledict[lines[1]+":"+lines[2]]=lines[0]
        return wegenefiledict

    '''
    /data4/compare/PMRA_make_full.txt
            affy_snp_id     chr_x   position_x      category        gene    rsid    chr_y   position_y      ref_allele      alt_allele
    0       Affx-2356540    10      100028647       GWAS grid               rs28366015      10      100028647       C       T
    1       Affx-2356580    10      100032749       GWAS grid               rs2862296       10      100032749       A       G
    '''
    def getPMRAfileDict(self,pmrafile):
        PMRADict={}
        with open(pmrafile,'r') as oper:
            for line in oper:
                lines=line.split('\t')
                PMRADict[lines[2]+":"+lines[3]]={"rsid":lines[6],'ref':lines[9]}
        return PMRADict

    def getVcfFileDict(self):
        '''
        :param vcffile: 输入vcf文件
        :return:key: chr:pos value:genotype 的字典
        '''
        vcfdict={}
        if os.path.isfile(self.vcffilename):
            vcf_reader = vcf.Reader(open(self.vcffilename, 'r'))
            sample=vcf_reader.samples[0]
            for record in vcf_reader:
                call=record.genotype(sample)
                vcfdict[record.CHROM[3:]+":"+str(record.POS)]={'ref':record.REF,'alt':record.ALT,'GT':record.genotype(sample)['GT'],'genotype':call.gt_bases}
            return vcfdict

        else:
            print("vcf文件不存在!")
            return vcfdict

    def getgenotype(self,gt,ref,alt):
        if gt=='0/1':
            return ref+alt
        elif(gt=='1/1'):
            return alt+alt
        elif(gt=='0/0'):
            return ref+ref
        else:
            return ref+ref

    def getVcfFileDictBySelf(self):
        '''
        #1 536134 . C T 0/1
        :param vcffile: 输入vcf文件
        :return:key: chr:pos value:genotype 的字典
        '''
        vcfdict={}
        if os.path.isfile(self.vcffilename):
            with open(self.vcffilename,'r') as vcfoper:
                for line in vcfoper:
                    if not line.startswith('#'):
                        lines=line.split(' ')
                        try:
                            vcfdict[lines[0]+":"+lines[1]]={'ref':lines[3],'GT':lines[5],'genotype':self.getgenotype(lines[5],lines[3],lines[4])}
                        except IndexError,e:
                            print(lines)
                            pass
            return vcfdict
        else:
            print("vcf文件不存在!")
            return vcfdict



    def printGen(self,vdict):
        if vdict['GT']=='0/0':
            vstr=vdict['genotype'].split('/')+['1','0','0']
        elif vdict['GT']=='0/1':
            vstr=vdict['genotype'].split('/')+['0','1','0']
        elif vdict['GT']=='1/1':
            vstr=vdict['genotype'].split('/')+['0','0','1']
        else:
            vstr=[vdict['ref']]*2+['1','0','0']
        return '\t'.join(vstr)

    def get138WithchromPos(self,chrom,pos):
        self.__cursor.execute("SELECT ref,rsid FROM T_DBSNP_HG19_138 WHERE CHROM='%s' and POS = '%s'"%(chrom,pos))
        results = self.__cursor.fetchone()
        if results!=None and results!="":
            return results[0],results[1]

    def get144WithchromPos(self,chrom,pos):
        self.__cursor.execute("SELECT ref,rsid FROM T_DBSNP_144_SEARCH WHERE CHROM='%s' and POS = '%s'"%(chrom,pos))
        results = self.__cursor.fetchone()
        if results!=None and results!="":
            return results[0],results[1]

    def get23andmeWithchromPos(self,chrom,pos):
        self.__cursor.execute("SELECT ref,rsid FROM t_23andme_addrefs WHERE chromosome='%s' and position = '%s'"%(chrom,pos))
        results = self.__cursor.fetchone()
        if results!=None and results!="":
            return results[0],results[1]

    def _23andme_left_join_vcf_gen(self,andmefile):
        '''
        6 rs499899 20001747 G G 1 0 0
        6 rs490655 20042503 G G 1 0 0
        6 rs558825 20046786 T T 1 0 0
        '''
        outfileoper=open('/tmp/23andme.15R4487.txt.TXT','w')
        if os.path.isfile(andmefile):
            andmedict=self.get23andmeFiledict(andmefile)
            vcfdict=self.getVcfFileDictBySelf()
            for k,v in andmedict.items():
                chrom,pos=k.split(':')
                if k in vcfdict:
                    outfileoper.writelines('\t'.join([chrom,v,pos,self.printGen(vcfdict[k])])+'\n')
                else:
                    try:
                        ref,rsid=self.get138WithchromPos(chrom,pos)
                        outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                    except Exception,e:
                        try:
                            ref,rsid=self.get144WithchromPos(chrom,pos)
                            outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                        except Exception,e:
                            try:
                                ref,rsid=self.get23andmeWithchromPos(chrom,pos)
                                outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                            except Exception,e:
                                print "ERROR:",chrom,pos
                                pass
        else:
            print("23andme 文件不存在")
        outfileoper.close()

    def genPrint(self):
        pass

    def wegene_left_join_vcf_gen(self,wegenefile):
        outfileoper=open('/tmp/wegene.15R4487.txt.TXT','w')
        if os.path.isfile(wegenefile):
            vvdict=self.getwegeneFiledict(wegenefile)
            vcfdict=self.getVcfFileDictBySelf()
            for k,v in vvdict.items():
                chrom,pos=k.split(':')
                if k in vcfdict:
                    outfileoper.writelines('\t'.join([chrom,v,pos,self.printGen(vcfdict[k])])+'\n')
                else:
                    try:
                        ref,rsid=self.get138WithchromPos(chrom,pos)
                        outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                    except Exception,e:
                        try:
                            ref,rsid=self.get144WithchromPos(chrom,pos)
                            outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                        except Exception,e:
                            try:
                                ref,rsid=self.get23andmeWithchromPos(chrom,pos)
                                outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                            except Exception,e:
                                print("ERROR:"+chrom+"\t"+pos)
                                pass
        else:
            print("wegenefile 文件不存在")
        outfileoper.close()

    def checkrsid(self,rsid,chrom,pos):
        if not rsid.startswith('rs'):
            try:
                ref,rsid=self.get138WithchromPos(chrom,pos)
                return rsid
            except Exception,e:
                try:
                    ref,rsid=self.get144WithchromPos(chrom,pos)
                    return rsid
                except Exception,e:
                    return "-"
        else:
            return rsid
    def PMRA_left_join_vcf_gen(self,PMRAfile):
        outfileoper=open('/tmp/PMRA.15R4487.txt.TXT','w')
        if os.path.isfile(PMRAfile):
            vcfdict=self.getVcfFileDictBySelf()
            andmedict=self.getPMRAfileDict(PMRAfile)
            for k,v in andmedict.items():
                chrom,pos=k.split(':')
                if k in vcfdict:
                    outfileoper.writelines('\t'.join([chrom,self.checkrsid(v[rsid],chrom,pos),pos,self.printGen(vcfdict[k])])+'\n')
                else:
                    ref=''
                    rsid=''
                    try:
                        ref,rsid=self.get138WithchromPos(chrom,pos)
                        outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                    except Exception,e:
                        try:
                            ref,rsid=self.get144WithchromPos(chrom,pos)
                            outfileoper.writelines('\t'.join([chrom,rsid,pos,ref,ref,'1','0','0'])+'\n')
                        except Exception,e:
                            print("ERROR:"+chrom+":"+pos)
                            pass
        else:
            print("23andme 文件不存在")
        outfileoper.close()



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

if __name__ == '__main__':
    # vcf_reader = vcf.Reader(open('/data4/compare/15R4487.raw.vcf.2', 'r'))
    # sample=vcf_reader.samples[0]
    # for record in vcf_reader:
    #     print(record.CHROM,record.POS,record.ID,record.REF,record.ALT,record.genotype(sample).gt_bases,record.genotype(sample)['GT'])
    k=VcfComvertToOtherFile("/data4/compare/15R4487.raw.vcf.short")
    #k.PMRA_left_join_vcf_gen('/data4/compare/PMRA_make_full.txt')
    #k._23andme_left_join_vcf_gen('/data4/compare/genome_wang_mian_Full_20151021231213.txt.bak')
    #k.wegene_left_join_vcf_gen('/data4/compare/wegeneRawData.txt')
    k.getSplitFile('/data4/share/impute/15R4487','/data4/compare/23andme_genome.gen')
