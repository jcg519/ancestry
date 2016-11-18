import re
import time
import MySQLdb
import optparse

class GenFile():
	def __init__(self,filename):
		self.filename = filename
		try:
			self.__db = MySQLdb.connect(host='192.168.30.252', port=3306, db='gendb', user='dna', passwd='dna',
										charset='utf8')
			self.__cursor = self.__db.cursor()
		except MySQLdb.Error, e:
			self.error_code = e.args[0]
			error_msg = 'MySQL error! ', e.args[0], e.args[1]
			print error_msg


if __name__ == "__main__":
	usage = "splitgen [ -i <genfule> -o <outdir>]"
	opter = optparse.OptionParser(usage)
	opter.add_option("-i", "--input",  dest="input", help=u"原始文件")
	opter.add_option("-o", "--output", dest="output", help=u"输出gen文件")
	opt, args = opter.parse_args()
	inputfile = opt.input
	outname = opt.output
	db=MySQLdb.connect(host='192.168.30.252',port=3306,db='gendb',user='dna',passwd='dna',charset='utf8')
#	cursor = db.cursor()
	a=open(inputfile,'r')
	outfileoper = open(outname, 'w')
	for line in a.readlines():
		if "#" not in line:
			res=re.split('\t',line)
			rsid=res[0]
			chrom=res[1]
			pos=res[2]
			genetype=res[3].strip()
			sql= "SELECT ref,rsid FROM T_DBSNP_144_SEARCH WHERE CHROM='%s' and POS = '%s' " %(chrom,pos)
			cursor.execute(sql)
			results = cursor.fetchall()
			for i in  results:
				ref=i[0]
				rs=i[1]
				if 'rs' in rsid:
					if genetype[0]==genetype[1]:
						if genetype[0] !='-':
#							chrom,rsid,pos,genetype[0],genetype[1],1,0,0
							outfileoper.writelines('\t'.join([chrom, rsid, pos, genetype[0],genetype[1], '1', '0', '0']) + '\n')
						else:
#							print chrom,rsid,pos,ref,ref,1,0,0
							outfileoper.writelines('\t'.join([chrom, rsid, pos, ref, ref, '1', '0', '0']) + '\n')
					else:
#						print chrom,rsid,pos,genetype[0],genetype[1],0,1,0
						outfileoper.writelines('\t'.join([chrom, rsid, pos, genetype[0], '0', '0', '1', '0']) + '\n')
				else:
					if genetype[0]==genetype[1]:
						if genetype[0] !='-':
#							print chrom,rs,pos,genetype[0],genetype[1],1,0,0
							outfileoper.writelines('\t'.join([chrom, rsid, pos, genetype[0], genetype[1], '1', '0', '0']) + '\n')
						else:
#							print chrom,rs,pos,ref,ref,1,0,0
							outfileoper.writelines('\t'.join([chrom, rsid, pos, ref, '0', '1', '0', '0']) + '\n')
					else:
#						print chrom,rs,pos,genetype[0],genetype[1],0,1,0
						outfileoper.writelines('\t'.join([chrom, rsid, pos, genetype[0], genetype[1], '0', '1', '0']) + '\n')

