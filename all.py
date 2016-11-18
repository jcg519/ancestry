__author__ = 'similarface'
import os
import sys
import optparse

class GetImpute:

    def getimpute(self,output):
        os.system(" cd "+ output)
        allsh=open('all.sh','w')
        for i in range(1,23):
            for j in os.listdir(output):
                if "chr"+str(i)+"_" in j:
                    shline=("/data4/share/program/impute_v2.3.2_x86_64_static/impute2 -h /data4/share/impute/ref/1000GP_Phase3_chr"+str(i)+
                  ".hap.gz  -l  /data4/share/impute/ref/1000GP_Phase3_chr"+str(i)+
                  ".legend.gz -m  /data4/share/impute/ref/genetic_map_chr"+str(i)+"_combined_b37.txt -g "+
                  os.path.join(output,j)+" -int "+
                  str(5000000*(int(j[j.index("_")+1:-4])-1))+" "+str(int(j[j.index("_")+1:-4])*5000000)+
                  " -Ne 20000  -o "+ output + str(j)+".imputed"+"\n")
                    allsh.write(shline)



if __name__=='__main__':
	usage = "splitgen [ -o <outputdir> "
	opter = optparse.OptionParser(usage)
	opter.add_option("-o", "--outputdir",dest="outputdir", help="the impute output dir")
	opt, args = opter.parse_args()
	outputdir = opt.outputdir
	impu=GetImpute()
	impu.getimpute(outputdir)    
