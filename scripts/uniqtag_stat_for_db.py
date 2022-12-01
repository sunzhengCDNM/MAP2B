#!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip
import threading
from datetime import datetime
from math import *
bindir = os.path.abspath(os.path.dirname(__file__))
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/11/20 02:31:28'
__version__ = '1.0.0'
############################################ main ##################################################
def report(level,info):
	date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	if level == "ERROR":
		sys.stderr.write("{0} - {1} - ERROR - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	elif level == "INFO":
		sys.stdout.write("{0} - {1} - INFO - {2}\n".format(date_now,os.path.basename(__file__),info))
	elif level == "DEBUG":
		sys.stdout.write("{0} - {1} - DEBUG - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	return()

def check_file(file):
	if os.path.exists(file):
		return(os.path.abspath(file))
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def check_dir(dir):
	dir = os.path.abspath(dir)
	if not os.path.exists(dir):
		os.system("mkdir -p {0}".format(dir))
		info = "mkdir: {0}".format(dir)
		report("INFO",info)
	return(dir)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-d',help='database file, CjePI.species.uniq.fa.gz',dest='database',type=str,required=True)
	parser.add_argument('-c',help='classify file, abfh_classify_with_speciename.txt.gz',dest='classify',type=str,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)
# reading classify file
	report('INFO', 'Start reading classify file')
	gcf_tax_dic = {} # {gcf:tax}
	with gzip.open(check_file(args.classify), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			gcf_tax_dic.setdefault(tmp[0], '\t'.join(tmp[1:8]))
	report('INFO', 'Complete reading classify file')
# reading database
	report('INFO', 'Start reading database')
	tag_gcf_dic = {} # {tag: [gcf]}
	gcf_theo_tag_num_dic = {} # {gcf: uniq_tag_num}
	tax_theo_tag_num_dic = {} # {tax: uniq_tag_num}
	db_file = check_file(args.database)
	enzyme = db_file.split('/')[-1].split('.')[0]
	with gzip.open(db_file, 'rt') as IN:
		for line in IN:
			if line.startswith('>'):
				line = line.strip()
				gcf = line.strip('>').split('|')[0]
			else:
				tag_gcf_dic.setdefault(line.strip(), []).append(gcf)
	for tag, gcf_list in tag_gcf_dic.items():
		for gcf in gcf_list:
			gcf_theo_tag_num_dic.setdefault(gcf, []).append(1)
	for gcf, t in gcf_theo_tag_num_dic.items():
		sum_t = sum(t)
		gcf_theo_tag_num_dic[gcf] = sum_t
		tax_theo_tag_num_dic.setdefault(gcf_tax_dic[gcf], []).append(sum_t)
	for tax, t in tax_theo_tag_num_dic.items():
		tax_theo_tag_num_dic[tax] = round(sum(t)/len(t), 4)
	with open(args.output, 'w') as OUT:
		for gcf, gcf_uniq_count in gcf_theo_tag_num_dic.items():
			OUT.write('{}\t{}\t{}\t{}\n'.format(gcf, gcf_tax_dic[gcf], gcf_theo_tag_num_dic[gcf], tax_theo_tag_num_dic[gcf_tax_dic[gcf]]))

if __name__=="__main__":
	main()
	info = "finish!"
#	report("INFO",info)
