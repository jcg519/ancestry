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
    parser.add_option("-m", "--name1", dest = "name1", help = "normal_name")
    parser.add_option("-n", "--name2", dest = "name2", help = "tumor_name")
    parser.add_option("-c", "--cancerType", dest = "cancer", default="lung",help = "lung cancer panel")
    parser.add_option("-r", "--ref", dest = "ref", default = "/data4/ctDNA_pepeline/ctDNA_ref/hg19.fa", help = "hg19 reference genome")
    parser.add_option("-t", "--thread", dest = "t", default = "20", help = "number of threads")
    parser.add_option("-o", "--output", dest = "output", default = "/data4/ctDNA_pepeline/data", help = "output file")
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

def annovar(snv_result, snv_cosmic):
	
	
	f = csv.DictReader(open("/data4/ctDNA_pepeline/ctDNA_ref/CosmicCompleteExport.tsv", "r"), delimiter="\t")
	cosmic = {}
	for i in f:
		if i["Gene name"] and i["Mutation CDS"]:
			gene_name=re.split("_",i["Gene name"])
			if re.search(">", i["Mutation CDS"]):
				tmp1 = re.search("(\d+)", i["Mutation CDS"])
				tmp2 = re.sub("(\d+)","", i["Mutation CDS"])
				tmp3 = re.sub(">", tmp1.groups()[0], tmp2)
				tmp = gene_name[0] +".*"+tmp3
			else:
				tmp = gene_name[0]+".*"+i["Mutation CDS"]
			if tmp not in cosmic:
				cosmic[tmp]=i["Mutation ID"]

	f_w = open(snv_cosmic, "a")
	header = ["gene","exon","acid","gDNA_varvant(%)" ,"ctDNA_varvant(%)","cosmic_id"]
	cosmic_write=csv.writer(f_w)
	cosmic_write.writerow(list(header))
	vcf = open(snv_result, "r")
	for i in vcf:
		if '#' not in i:
			vcflist=re.split("\t", i)
			normal = re.match("(.*?):.:.*?:(.*?):(.*?):(.*?)%:.*?,.*?,.*?,.*?", vcflist[9])  #####  normal sample
			tumor = re.match("(.*?):.:.*?:(.*?):(.*?):(.*?)%:.*?,.*?,.*?,.*?", vcflist[10])  #####  cancer sample
			acid = re.findall("AAChange.refGene=(.*?);",vcflist[7])
			Band = re.findall("cytoBand=(.*?);",vcflist[7])
			Func = re.findall("Func.refGene=(.*?);",vcflist[7])
			j = re.split(",", acid[0])
			ctdna = tumor.groups()[3]
			gdna  = normal.groups()[3]

			list_ll=[vcflist[0],vcflist[1],vcflist[3],vcflist[4],Band[0],Func[0],j[0],normal.groups()[0],normal.groups()[1],normal.groups()[2],tumor.groups()[0],tumor.groups()[1],tumor.groups()[2],gdna,ctdna]
			if list_ll[5]=="exonic":
				if float(list_ll[14]>0.05):
					if list_ll[6] !=".":
						res = re.match("^(.*?):.*?:(.*?):(c.*?):(.*?)$", list_ll[6])
						if res.groups()[0]+".*"+res.groups()[2] in cosmic:
							idd=cosmic[res.groups()[0]+".*"+res.groups()[2]]
					
							filerow=[str(res.groups()[0]), str(res.groups()[1]+':'+res.groups()[2]), str(res.groups()[3]), str(list_ll[13]),str(list_ll[14]), str(str("COSM")+''+idd)]
					
#							print filerow
							cosmic_write.writerow(filerow)
						else:
							idd=' '
							filerow=[str(res.groups()[0]), str(res.groups()[1]+':'+res.groups()[2]), str(res.groups()[3]), str(list_ll[13]), str(list_ll[14]), str(idd)]

#							print filerow
							cosmic_write.writerow(filerow)
        f_w.close()

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
#    if not os.path.exists(sam1):
#        logging.info("start bwa mem")
#        os.system("bwa mem -t " + t + " " + ref + " " + input1 +  " " + input2 + " > " + sam1)
    if not os.path.exists(bam_sort1):
	logging.info("start samtools sort")
	os.system("bwa mem -t"+ t +" " + ref + " " + input1 + " " + input2 + " |" + "samtools view -Sh -" + " | " + "samtools view -uS -" + " |" + "samtools sort -m 2000000000 -  "+ bam_sort1)
    logging.info("============================")
    if not os.path.exists(bam_mark1):
	logging.info("start Picard duplicate mark and remove")
	os.system("/data4/ctDNA_pepeline/normal_tumor/dedup.py -1 " + bam_sort1 + " -o " + bam_mark1 )
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
#    if not os.path.exists(sam2):
#        logging.info("start bwa mem")
#        os.system("bwa mem -t " + t + " " + ref + " " + input1 +  " " + input2 + " > " + sam2)
#    logging.info("============================")
    if not os.path.exists(bam_sort2):
	logging.info("start samtools sort")
	os.system("bwa mem -t"+ t +" " + ref + " " + input1 + " " + input2 + " |" + "samtools view -Sh -" + " | " + "samtools view -uS -" + " |" + "samtools sort -m 200
    0000000 -  " + bam_sort2)
    logging.info("============================")
    if not os.path.exists(bam_mark2):
	logging.info("start Picard duplicate mark and remove")
	os.system("/data4/ctDNA_pepeline/normal_tumor/dedup.py -1 " + bam_sort2 + " -o " + bam_mark2 )
    logging.info("============================") 
    if not os.path.exists(bam_rg2):
	logging.info("start Picard add RG")
	os.system("java -Xmx10g -jar /home/cjiang/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=" + bam_mark2 + " OUTPUT=" + bam_rg2 + " RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample003")
    logging.info("============================")
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
    annovar_snp_vcf_prefix   = output + "/" + name + "_snp_annovar"
    annovar_indel_vcf_prefix   = output + "/" + name + "_indel_annovar"
    annovar_snp_vcf   = output + "/" + name + "_snp_annovar.hg19_multianno.vcf"
    annovar_indel_vcf   = output + "/" + name + "_indel_annovar.hg19_multianno.vcf"
    snp_new_vcf   = output + "/" + name + "_snp_annovar_new.hg19_multianno.vcf"
    indel_new_vcf   = output + "/" + name + "_indel_annovar_new.hg19_multianno.vcf"
    snp_result_vcf   = output + "/" + name + "_snp_annovar_result.hg19_multianno.vcf"
    indel_result_vcf   = output + "/" + name + "_indel_annovar_result.hg19_multianno.vcf"
    snp_cosmic         = output + "/" + name + "_snp_cosmic.vcf"
    indel_cosmic       = output + "/" + name + "_indel_cosmic.vcf"
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',datefmt = '%a, %d %b %Y %H:%M:%S',filename = log,filemode = 'w')
    if not os.path.exists(snp):
	    logging.info("start compare normal tumor")
	    os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar somatic ' + normal + ' ' + tumor +  ' ' + normal_tumor + " --min-var-freq 0.0005 --tumor-purity 0.001 --output-vcf 1")
    logging.info("===========================")
    if not os.path.exists(snp_filter):
		logging.info("start somaticFilter")
		os.system('java -Xmx10g -jar /home/cjiang/tools/varscan/VarScan.v2.4.0.jar somaticFilter '+ snp + " --indel-file " + indel + " --min-var-freq 0.000001 --p-value 1 --output-file " + snp_filter)
    logging.info("===========================")
    if not os.path.exists(annovar_snp_vcf):
        logging.info("start get annotation")
        os.system('/data1/annovar/table_annovar.pl ' + snp_filter + ' ' + '/data1/annovar/humandb/ -buildver hg19 -out ' + annovar_snp_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
        os.system('/data1/annovar/table_annovar.pl ' + indel + ' ' + '/data1/annovar/humandb/ -buildver hg19 -out ' + annovar_indel_vcf_prefix  + ' ' + '-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput')
    if not os.path.exists(snp_new_vcf):
	    logging.info("start remove dbSNP")
	    os.system("/data4/ctDNA_pepeline/normal_tumor/vcf_filter.py -1 " + annovar_snp_vcf + " -g 0.001 -o " + snp_new_vcf)
	    os.system("/data4/ctDNA_pepeline/normal_tumor/vcf_filter.py -1 " + annovar_indel_vcf + " -g 0.001 -o " + indel_new_vcf)
    logging.info("===========================")
    if not os.path.exists(snp_result_vcf):
	    logging.info("start extract result")
	    gene(snp_new_vcf, snp_result_vcf, cancer)
	    gene(indel_new_vcf, indel_result_vcf, cancer)
    logging.info("===========================")
    if not os.path.exists(snp_cosmic):
        logging.info("start extract result")
        annovar(snp_result_vcf, snp_cosmic)
        annovar(indel_result_vcf, indel_cosmic)
    logging.info("===========================")


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
    t1 = threading.Thread( target=getBam3, args=(options.input1, options.input2, options.name1, options.output, options.t, options.ref) )
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
