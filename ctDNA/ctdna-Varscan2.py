#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os,sys
import logging
import commands
from optparse import OptionParser
import time
import gzip
import csv
import re


def parseCommand():
	usage = "usage: %prog <-1 input1> <-2 input2>  <-n name> <-o outputDir>"
	version = "%prog 1.0"
	parser = OptionParser(usage = usage, version = version)
	parser.add_option("-1", "--input1", dest = "input1", help = "the read1 file")
	parser.add_option("-2", "--input2", dest = "input2", help = "the read2 file")
	parser.add_option("-n", "--name", dest = "name", help = "输出文件的前缀名".decode("utf8"))
	parser.add_option("-c", "--cancerType", dest = "cancer", default="lung",help = "肿瘤类型 ".decode("utf8"))
	parser.add_option("-r", "--ref", dest = "ref", default = "/data4/ctDNA_pepeline/ctDNA_ref/hg19.fa", help = "hg19 reference genome")
	parser.add_option("-t", "--thread", dest = "t", default = "20", help = "number of threads")
	parser.add_option("-o", "--output", dest = "output", default = "/data4/ctDNA_pepeline/data", help = "输出目录".decode("utf8"))
	return parser.parse_args()

def gene(snv, snv_result, cancer='lung'):
	genes = []
	csvfile = open('/data4/ctDNA_pepeline/ctDNA_ref/panel.csv', "r")
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

if __name__ == "__main__":
	time1 = time.time()
	(options, args) = parseCommand()
	if options.input1 == None:
		print ("input1(reads1) is required, use -1 <file name> to specify")
		print ("see -h for help")
		sys.exit(-1)
	if options.input2 == None:
		print ("input2(reads2) is required, use -2 <file name> to specify")
		print ("see -h for help")
		sys.exit(-1)
	if options.name == None:
		print ("output name is required, use -n <file prefix name> to specify")
		print ("see -h for help")
		sys.exit(-1)
	sam            = options.output + "/" + options.name + ".sam"
	bam_sort       = options.output + "/" + options.name + "_sort.bam"
	bam_mark       = options.output + "/" + options.name + "_mark.bam"
	bam_mark_txt   = options.output + "/" + options.name + "_metrics.txt"
	bam_rg         = options.output + "/" + options.name + "_rg.bam"
	bam_rg_index   = options.output + "/" + options.name + "_rg.bai"
	bam_rg_new     = options.output + "/" + options.name + "_bam_rg_new.bam"
	mpileup        = options.output + "/" + options.name + ".mpileup"
	snp            = options.output + "/" + options.name + "_varscan.snp.vcf"
	indel          = options.output + "/" + options.name + "_varScan.indel.vcf"
	snp_filter     = options.output + "/" + options.name + "_snp.filter.vcf"
	indel_filter   = options.output + "/" + options.name + "_indel.filter.vcf"
	annovar_snp_vcf_prefix    = options.output + "/" + options.name + "_snp_annovar"
	annovar_indel_vcf_prefix  = options.output + "/" + options.name + "_indel_annovar"
	annovar_snp_vcf           = options.output + "/" + options.name + "_snp_annovar.hg19_multianno.vcf"
	annovar_indel_vcf         = options.output + "/" + options.name + "_indel_annovar.hg19_multianno.vcf"
	annovar_new_snp_vcf       = options.output + "/" + options.name + "_annovar_new_snp.hg19_multianno.vcf"
	annovar_new_indel_vcf     = options.output + "/" + options.name + "_annovar_new_indel.hg19_multianno.vcf"
	result_snp_vcf     = options.output + "/" + options.name + "_annovar_snp_result.vcf"
	result_indel_vcf   = options.output + "/" + options.name + "_annovar_indel_result.vcf"
	log = options.output + "/" + options.name + ".log"
	logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
	if not os.path.exists(sam):
	    logging.info("start bwa mem")
	    os.system("bwa mem -t " + options.t + " " + options.ref + " " + options.input1 +  " " + options.input2 + " > " + sam)
	logging.info("============================")
	if not os.path.exists(bam_sort):
	    logging.info("start Picard sort")
	    os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/SortSam.jar INPUT=" + sam + " OUTPUT=" + bam_sort + " SORT_ORDER=coordinate")
	logging.info("============================")
	if not os.path.exists(bam_mark):
	    logging.info("start Picard duplicate mark")
	    os.system("/data4/ctDNA_pepeline/normal_tumor/dedup.py -1 " + bam_sort + " -o " + bam_mark)
	logging.info("============================")
	if not os.path.exists(bam_rg):
		logging.info("start Picard add RG")
		os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=" + bam_mark + " OUTPUT=" + bam_rg + " RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 "+" RGSM="+options.name)
	logging.info("============================")
	if not os.path.exists(bam_rg_index):
	    logging.info("start Picard build index")
	    os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/BuildBamIndex.jar INPUT=" + bam_rg)
	logging.info("============================")
	if not os.path.exists(mpileup):
	    logging.info("start generate pileup")
	    os.system("samtools view -b -u -q 1 " + bam_rg + " > " + bam_rg_new)
	    os.system("samtools mpileup -f " + options.ref + " " + bam_rg_new + "|" + "awk '{if($4 != 0) print $0}'"+ " > " + mpileup)
	logging.info("============================")
	if not os.path.exists(snp):
		logging.info("start varscan2")
		os.system('java -Xmx40g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar mpileup2snp' + ' '+ mpileup + ' '+ '--min-coverage 5 --min-var-freq 0.01 --p-value 0.05  -output-vcf 1' + ">" + snp)
		os.system('java -Xmx40g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar mpileup2indel' + ' ' + mpileup  + ' '+ '--min-coverage 5 --min-var-freq 0.05 --p-value 0.10  -output-vcf 1' + ">" + indel)
	logging.info("============================")
	if not os.path.exists(snp_filter):
		logging.info("start filter")
		os.system('java -Xmx40g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar filter' + ' '+  snp  + ' '+ '--indel-file'+ ' ' + indel + ' '+ '--output-file' +' '+ snp_filter)
		os.system('java -Xmx40g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar filter' + ' '+  indel  +' '+ '--min-reads 1 --min-var-freq 0.15 --p-value 0.05 --output-file' + ' ' +  indel_filter)
	logging.info("============================")
	if not os.path.exists(annovar_snp_vcf):
		logging.info("start annovar")
		os.system("/data1/Shark1.1/Tools/annovar/table_annovar.pl" + ' '+ snp_filter  +' '+ '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out' + ' '+  annovar_snp_vcf_prefix + " -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput")
		os.system("/data1/Shark1.1/Tools/annovar/table_annovar.pl" + ' '+ indel_filter +' '+ '/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out' +' ' + annovar_indel_vcf_prefix + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
	logging.info("============================")
	if not os.path.exists(annovar_new_snp_vcf):
		logging.info("start annovar")
		os.system("/data4/ctDNA_pepeline/normal_tumor/vcf_filter.py -1 " + ' '+ annovar_snp_vcf + " -g 0.001 -o " + ' ' + annovar_new_snp_vcf)
		os.system("/data4/ctDNA_pepeline/normal_tumor/vcf_filter.py -1 " + ' '+ annovar_indel_vcf + " -g 0.001 -o " + ' ' + annovar_new_indel_vcf)
	logging.info("============================")
	if not os.path.exists(result_indel_vcf):
		logging.info("start filter by panel gene")
		gene(annovar_new_snp_vcf, result_snp_vcf, options.cancer)
		gene(annovar_new_indel_vcf, result_indel_vcf, options.cancer)
	logging.info("============================")

	time2 = time.time()
	print options.ref
#	print options.t == "20"
	print 'Time used: ' + str(time2-time1)
