#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os,sys
import logging
from optparse import OptionParser
import time
import gzip
import csv
import re
import threading


def parseCommand():
    usage = "usage: %prog <-1 input1> <-2 input2>  <-n name> <-o outputDir>"
    version = "%prog 1.0"
    parser = OptionParser(usage = usage, version = version)
    parser.add_option("-1", "--input1", dest = "input1", help = "the read1 file")
    parser.add_option("-2", "--input2", dest = "input2", help = "the read2 file")
    parser.add_option("-3", "--input3", dest = "input3", help = "the read3 file")
    parser.add_option("-4", "--input4", dest = "input4", help = "the read4 file")
    parser.add_option("-m", "--name1", dest = "name1", help = "normal输出文件的前缀名".decode("utf8"))
    parser.add_option("-n", "--name2", dest = "name2", help = "tumor输出文件的前缀名".decode("utf8"))
    parser.add_option("-c", "--cancerType", dest = "cancer", default="lung",help = "肿瘤类型,比如肺癌lung 乳腺癌breast 直肠癌rectum 胃癌stomach 肝癌liver 以及48genes".decode("utf8"))
    parser.add_option("-r", "--ref", dest = "ref", default = "/data3/20160125WholeGenomeAnalysis/Data/Reference/hg19.fa", help = "hg19 reference genome")
    parser.add_option("-t", "--thread", dest = "t", default = "20", help = "number of threads")
    parser.add_option("-o", "--output", dest = "output", default = "/data3/cjiang/data", help = "输出目录".decode("utf8"))
    return parser.parse_args()

def gene(snv, snv_result, cancer='lung'):
    genes = []
    csvfile = open('/home/cjiang/panel.csv', "r")
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row[cancer]:
            genes.append(row[cancer])
    csvfile.close()
    fd_snv    = open(snv, "r")
    fd_result = open(snv_result, "a")
    for line in fd_snv.readlines():
        for gene in genes:
            if re.search("[=,]"+gene+"[,;]", line):
                fd_result.write(line)
    fd_snv.close()
    fd_result.close()
    return 1 


def getBam3(input1, input2, name1, output, t, ref):
    sam1            = output + "/" + name1 + ".sam"
    bam_sort1       = output + "/" + name1 + "_sort.bam"
    bam_mark1       = output + "/" + name1 + "_mark.bam"
    bam_mark_index1 = output + "/" + name1 + "_mark.bai"
    bam_rg_new1     = output + "/" + name1 + "_rg_new.bam"
    pileup1         = output + "/" + name1 + ".pileup"
    log = output + "/" + name1 + ".log"
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
    if not os.path.exists(sam1):
        logging.info("start bwa mem")
        os.system("bwa mem -t " + t + " " + ref + " " + input1 +  " " + input2 + " > " + sam1)
    if not os.path.exists(bam_sort1):
	logging.info("start samtools sort")
	os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/SortSam.jar INPUT=" + sam1 + " OUTPUT=" + bam_sort1 + " SORT_ORDER=coordinate")
    logging.info("============================")
    if not os.path.exists(bam_mark1):
	logging.info("start Picard duplicate mark and remove")
#	os.system("java -Xmx10g -jar /tools/GATK/picard/MarkDuplicates.jar M=" + bam_mark_txt + " I=" + bam_sort + " O=" + bam_mark + " REMOVE_DUPLICATES=true")
	os.system("/home/cjiang/pepeline/dedup.py -1 " + bam_sort1 + " -o " + bam_mark1 )
    logging.info("============================")
    if not os.path.exists(bam_mark_index1):
	logging.info("start Picard build index")
        os.system("samtools index " + bam_mark1)
    logging.info("============================")
    if not os.path.exists(pileup1):
	logging.info("start generate pileup")
	os.system("samtools view -b -u -q 1 " + bam_mark1 + " > " + bam_rg_new1)
	os.system("samtools mpileup -f " + ref + " " + bam_rg_new1 + " > " + pileup1)
    logging.info("============================")
    return 1

def getBam2(input1, input2, name2, output, t, ref):
    sam2                   = output + "/" + name2 + ".sam"
    bam2                   = output + "/" + name2 + ".bam"
    bam_sort2              = output + "/" + name2 + "_sort.bam"
    bam_mark2              = output + "/" + name2 + "_mark.bam"
    bam_mark_index2        = output + "/" + name2 + "_mark.bai"
    bam_rg_new2            = output + "/" + name2 + "_rg_new.bam"
    bam_mark_txt2          = output + "/" + name2 + "_metrics.txt"
    bam_rg2                = output + "/" + name2 + "_rg.bam"
    bam_rg_index2          = output + "/" + name2 + "_rg.bai"
    bam_rg_new2            = output + "/" + name2 + "_rg_new.bam"
    GATK_realn_intervals2  = output + "/" + name2 + "_Realn.intervals"
    bam_realn2             = output + "/" + name2 + "_Realn.bam"
    bam_baserecal_table2   = output + "/" + name2 + "_baserecal.table"
    bam_baserecal2         = output + "/" + name2 + "_baserecal.bam"
    pileup2                = output + "/" + name2 + ".pileup"
    log = output + "/" + name2 + ".log"
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
    if not os.path.exists(sam2):
        logging.info("start bwa mem")
        os.system("bwa mem -t " + t + " " + ref + " " + input1 +  " " + input2 + " > " + sam2)
    logging.info("============================")
    if not os.path.exists(bam_sort2):
	logging.info("start Picard sort")
	os.system("java -Xmx10g -jar /home/cjiang/tools/picard-tools-1.119/SortSam.jar INPUT=" + sam2 + " OUTPUT=" + bam_sort2 + " SORT_ORDER=coordinate")
    logging.info("============================")
    if not os.path.exists(bam_mark2):
	logging.info("start Picard duplicate mark and remove")
	os.system("/home/cjiang/pepeline//dedup.py -1 " + bam_sort2 + " -o " + bam_mark2 )
    logging.info("============================") 
    if not os.path.exists(bam_rg2):
	logging.info("start Picard add RG")
	os.system("java -Xmx10g -jar /home/cjiang/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=" + bam_mark2 + " OUTPUT=" + bam_rg2 + " RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample003")
    logging.info("============================")
#    if not os.path.exists(bam_rg_index):
    if not os.path.exists(bam_mark_index2):
	logging.info("start Picard build index")
	os.system("java -Xmx10g -jar /home/cjiang/tools/picard-tools-1.119/BuildBamIndex.jar INPUT=" + bam_rg2)
    logging.info("============================")
    if not os.path.exists(pileup2):
	logging.info("start generate pileup")
	os.system("samtools view -b -u -q 1 " + bam_rg2 + " > " + bam_rg_new2)
	os.system("samtools mpileup -f " + ref + " " + bam_rg_new2 + " > " + pileup2)
    logging.info("============================")
    return 1


def normal_tumor_compare3(output, name, normal, tumor, cancer):
	normal_tumor = output + "/" + name + " "
	snp_filter   = output + "/" + name + "_varscan2_filter.snp.vcf"
	log        = output + "/" + name  + ".log"
	snp        = output + "/" + name + ".snp.vcf"
	indel      = output + "/" + name + ".indel.vcf"
	snp_hc     = output + "/" + name + ".snp.Somatic.hc.vcf"
	indel_hc   = output + "/" + name + ".indel.Somatic.hc.vcf"
	snp_loh_hc	  = output + "/" + name  + ".snp.LOH.hc.vcf"
	indel_loh_hc    = output + "/" + name  + ".indel.LOH.hc.vcf"
	germline_snp_hc =output + "/" + name  + ".snp.Germline.hc.vcf"	
	germline_indel_hc =output + "/" + name + ".indel.Germline.hc.vcf" 
	annovar_snp_vcf_prefix   = output + "/" + name + "_snp_annovar"
	annovar_indel_vcf_prefix   = output + "/" + name + "_indel_annovar"
	annovar_snp_loh_vcf_prefix   = output + "/" + name + "_snp_loh_annovar"
	annovar_indel_loh_vcf_prefix = output + "/" + name + "_indel_loh_annovar"
	annovar_snp_germline_vcf_prefix = output + "/" + name + "_snp_germline_annovar"
	annovar_indel_germline_vcf_prefix= output + "/" + name + "_indel_germline_annovar"
	annovar_snp_vcf   = output + "/" + name + "_snp_annovar.hg19_multianno.vcf"
	annovar_indel_vcf   = output + "/" + name + "_indel_annovar.hg19_multianno.vcf"
	annovar_snp_loh_vcf   = output + "/" + name + "_snp_loh_annovar.hg19_multianno.vcf"
	annovar_indel_loh_vcf   = output + "/" + name + "_indel_loh_annovar.hg19_multianno.vcf"
	annovar_snp_germline_vcf = output + "/" + name + "_snp_germline_annovar.hg19_multianno.vcf"
	annovar_indel_germline_vcf= output + "/" + name + "_indel_germline_annovar.hg19_multianno.vcf"
	snp_new_vcf   = output + "/" + name + "_snp_annovar_new.hg19_multianno.vcf"
	snp_loh_new_vcf   = output + "/" + name + "_snp_loh_annovar_new.hg19_multianno.vcf"
	indel_new_vcf   = output + "/" + name + "_indel_annovar_new.hg19_multianno.vcf"
	snp_new_germline_vcf = output + "/" + name + "_snp_germline_annovar_new.hg19_multianno.vcf"
	indel_loh_new_vcf   = output + "/" + name + "_indel_loh_annovar_new.hg19_multianno.vcf"
	indel_new_germline_vcf = output + "/" + name + "_indel_germline_annovar_new.hg19_multianno.vcf"
	snp_result_vcf   = output + "/" + name + "_snp_annovar_result.hg19_multianno.vcf"
	snp_loh_result_vcf  = output + "/" + name + "_snp_loh_annovar_result.hg19_multianno.vcf"
	indel_result_vcf   = output + "/" + name + "_indel_annovar_result.hg19_multianno.vcf"
	indel_loh_result_vcf   = output + "/" + name + "_indel_loh_annovar_result.hg19_multianno.vcf"
	snp_germline_result_vcf =output + "/" +name + "_snp_germline_result.hg19_multianno.vcf"
	indel_germline_result_vcf =output + "/" +name + "_indel_germline_result.hg19_multianno.vcf"
	logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
	if not os.path.exists(snp):
		logging.info("start compare normal tumor")
		os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar somatic ' + normal + ' ' + tumor +  ' ' + normal_tumor + " --min-var-freq 0.0005 --tumor-purity 0.001 --output-vcf 1")
	logging.info("===========================")
	if not os.path.exists(snp_hc):
		logging.info("start split normal tumor")
		os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar processSomatic' + ' ' + snp)
		os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar processSomatic' + ' ' + indel)
	logging.info("===========================")
	if not os.path.exists(snp_filter):
		logging.info("start somaticFilter")
		os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar somaticFilter '+ snp_hc + " --indel-file " + indel + " --min-var-freq 0.000001 --p-value 1 --output-file " + snp_filter)
	logging.info("===========================")
	if not os.path.exists(annovar_indel_loh_vcf):
		logging.info("start get annotation")
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + snp_filter + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_snp_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + indel_hc + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_indel_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + snp_loh_hc + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_snp_loh_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + germline_snp_hc + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_snp_germline_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + indel_loh_hc + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_indel_loh_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		os.system('/data1/Shark1.1/Tools/annovar/table_annovar.pl ' + germline_indel_hc + ' ' + '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out ' + annovar_indel_germline_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
		

	logging.info("===========================")
	if not os.path.exists(snp_new_vcf):

		logging.info("start remove dbSNP")
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_snp_vcf + " -g 0.001 -o " + snp_new_vcf)
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_snp_loh_vcf + " -g 0.001 -o " + snp_loh_new_vcf)
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_indel_vcf + " -g 0.001 -o " + indel_new_vcf)
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_indel_loh_vcf + " -g 0.001 -o " + indel_loh_new_vcf)
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_snp_germline_vcf + " -g 0.001 -o " + snp_new_germline_vcf)
		os.system("/home/cjiang/pepeline/vcf_filter.py -1 " + annovar_indel_germline_vcf + " -g 0.001 -o " + indel_new_germline_vcf)
	logging.info("===========================")
	if not os.path.exists(indel_germline_result_vcf):
		logging.info("start extract result")
		gene(snp_new_vcf, snp_result_vcf, cancer)
		gene(indel_new_vcf, indel_result_vcf, cancer)
		gene(snp_loh_new_vcf, snp_loh_result_vcf, cancer)
		gene(indel_loh_new_vcf, indel_loh_result_vcf, cancer)
		gene(snp_new_germline_vcf, snp_germline_result_vcf, cancer)
		gene(indel_new_germline_vcf, indel_germline_result_vcf, cancer)


if __name__ == "__main__":
    time1 = time.time()
    (options, args) = parseCommand()
    if options.input1 == None:
	print "input1(normal reads1) is required, use -1 <file name> to specify"
	print "see -h for help"
	sys.exit(-1)
    if options.input2 == None:
	print "input2(normal reads2) is required, use -2 <file name> to specify"
	print "see -h for help"
	sys.exit(-1)
    if options.input3 == None:
	print "input3(tumor reads3) is required, use -3 <file name> to specify"
	print "see -h for help"
	sys.exit(-1)
    if options.input4 == None:
	print "input4(tumor reads4) is required, use -4 <file name> to specify"
	print "see -h for help"
	sys.exit(-1)
    if options.name1 == None:
	print "output name1 is required, use -m <file prefix name> to specify"
	print "see -h for help"
	sys.exit(-1)
    if options.name2 == None:
	print "output name2 is required, use -n <file prefix name> to specify"
	print "see -h for help"
	sys.exit(-1)
    threads = []
    t1 = threading.Thread( target=getBam2, args=(options.input1, options.input2, options.name1, options.output, options.t, options.ref) )
    threads.append(t1)
    t2 = threading.Thread( target=getBam2, args=(options.input3, options.input4, options.name2, options.output, options.t, options.ref) )
    threads.append(t2)
    for t in threads:
        t.setDaemon(True)
	t.start()
    for t in threads:
	t.join()
    normal = options.output + "/" + options.name1 + ".pileup"
    tumor  = options.output + "/" + options.name2 + ".pileup"
    normal_tumor_compare3(options.output, options.name2, normal, tumor, options.cancer)
    time2 = time.time()
    print options.ref
    print 'Time used:' + str(time2-time1)
