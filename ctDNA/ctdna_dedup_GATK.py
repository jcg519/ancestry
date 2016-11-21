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
import multiprocessing


def parseCommand():
	usage = "usage: %prog <-1 input1> <-2 input2> <-n name> <-o outputDir>"
	version = "%prog 1.0"
	parser = OptionParser(usage = usage, version = version)
	parser.add_option("-1", "--input1", dest = "input1", help = "the read1 file")
	parser.add_option("-2", "--input2", dest = "input2", help = "the read2 file")
	parser.add_option("-n", "--name", dest = "name", help ="ouput name ")
	parser.add_option("-c", "--cancerType", dest = "cancer",help = "lung,breast ")
	parser.add_option("-r", "--ref", dest = "ref", default = "/data4/ctDNA_pepeline/ctDNA_ref/hg19.fa", help = "hg19 reference genome")
	parser.add_option("-t", "--thread", dest = "t", default = "20", help = "number of threads")
	parser.add_option("-o", "--output", dest = "output", default = "/data4/result", help = "输出目录".decode("utf8"))
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
	sam                  = options.output + "/" + options.name + ".sam"
	bam_sort             = options.output + "/" + options.name + "_sort.bam"
	bam_mark             = options.output + "/" + options.name + "_mark.bam"
	bam_mark_txt         = options.output + "/" + options.name + "_metrics.txt"
	bam_rg               = options.output + "/" + options.name + "_rg.bam"
	bam_rg_index         = options.output + "/" + options.name + "_rg.bai"
	GATK_realn_intervals = options.output + "/" + options.name + "_Realn.intervals"
	bam_realn            = options.output + "/" + options.name + "_Realn.bam"
	bam_baserecal_table  = options.output + "/" + options.name + "_baserecal.table"
	bam_baserecal        = options.output + "/" + options.name + "_baserecal.bam"
	vcf                  = options.output + "/" + options.name + "_HaplotyperCalled.vcf"
	annovar_vcf_prefix   = options.output + "/" + options.name + "_annovar"
	annovar_vcf          = options.output + "/" + options.name + "_annovar.hg19_multianno.vcf"
	annovar_new_vcf      = options.output + "/" + options.name + "_annovar_new.hg19_multianno.vcf"
	result_vcf           = options.output + "/" + options.name + "_annovar_result.hg19_multianno.vcf"
	log = options.output + "/" + options.name + ".log"
	logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
	if not os.path.exists(clean):
		logging.info("start qc")
		os.system("")
		logging.info("============================")

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
#		os.system("java -Xmx40g -jar /tools/GATK/picard/MarkDuplicates.jar M=" + bam_mark_txt + " I=" + bam_sort + " O=" + bam_mark)
        os.system("/data4/ctDNA_pepeline/normal_tumor/dedup.py -1 " + bam_sort + " -o " + bam_mark )
	logging.info("============================")
        if not os.path.exists(bam_rg):
		logging.info("start Picard add RG")
		os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=" + bam_mark + " OUTPUT=" + bam_rg + " RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 "+" RGSM="+options.name)
	logging.info("============================")
	if not os.path.exists(bam_rg_index):
		logging.info("start Picard build index")
		os.system("java -Xmx40g -jar /home/cjiang/tools/picard-tools-1.119/BuildBamIndex.jar INPUT=" + bam_rg)
	logging.info("============================")
	if not os.path.exists(GATK_realn_intervals):
		logging.info("start GATK Indel Realignment")
		os.system("/usr/bin/java -Xmx40g -jar /data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + options.ref + " -nt 30 -I " + bam_rg + " -o " + GATK_realn_intervals + " -known /data3/cjiang/ref/Mills_and_1000G_gold_standard.indels.hg19.vcf")
	logging.info("============================")
	if not os.path.exists(bam_realn):
		logging.info("start GATK Indel Realignment 2")
		os.system("/usr/bin/java -Xmx40g -jar /data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R " + options.ref + " -I " + bam_rg + " -targetIntervals " + GATK_realn_intervals + " -known /data3/cjiang/ref/Mills_and_1000G_gold_standard.indels.hg19.vcf -o " + bam_realn)
	logging.info("============================")
	if not os.path.exists(bam_baserecal_table):
		logging.info("start GATK BaseRecalibration")
		os.system("/usr/bin/java -Xmx40g -jar /data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 40 -R " + options.ref + " -knownSites /data3/cjiang/ref/dbsnp144hg19.vcf " +" -I " + bam_realn + " -o " + bam_baserecal_table)
	logging.info("============================")
	if not os.path.exists(bam_baserecal):
		logging.info("start GATK BaseRecalibration 2")
		os.system("/usr/bin/java -Xmx40g -jar /data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar -T PrintReads -nct 40 -R " + options.ref + " -I " + bam_realn + " -BQSR " + bam_baserecal_table + " -o " + bam_baserecal)
	logging.info("============================")
	if not os.path.exists(vcf):
		logging.info("start GATK HaplotypeCaller")
		os.system("/usr/bin/java -Xmx40g -jar /data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + options.ref + " -I " + bam_baserecal +  " -o " + vcf + " -stand_call_conf 30 -stand_emit_conf 10 -minPruning 3")
	logging.info("============================")
	if not os.path.exists(annovar_vcf):
		logging.info("start annovar")
		os.system("/data1/Shark1.1/Tools/annovar/table_annovar.pl " + ' ' +vcf + ' ' + "/data1/Shark1.1/Tools/annovar/humandb/ -buildver hg19 -out " +  annovar_vcf_prefix + " -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput")
	logging.info("============================")
	if not os.path.exists(annovar_new_vcf):
		logging.info("start annovar")
		os.system("/data4/ctDNA_pepeline/normal_tumor/vcf_filter.py -1 " + annovar_vcf + " -g 0.001 -o " + annovar_new_vcf)
	logging.info("============================")
	if not os.path.exists(result_vcf):
		logging.info("start filter by panel gene")
		gene(annovar_new_vcf, result_vcf, options.cancer)
	logging.info("============================")

	time2 = time.time()
	print options.ref
	print 'Time used: ' + str(time2-time1)
