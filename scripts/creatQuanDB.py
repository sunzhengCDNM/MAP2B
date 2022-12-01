#!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip
from datetime import datetime
from itertools import chain
bindir = os.path.abspath(os.path.dirname(__file__))
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/11/24 22:10:41'
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
	parser.add_argument('-d',help='db file',dest='database',type=str,required=True)
	parser.add_argument('-c',help='classfy file',dest='classfy',type=str,required=True)
	parser.add_argument('-p',help='pred file',dest='pred',type=str,required=True)
	parser.add_argument('-o',help='new db file',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)
	spe_list = []
	with open(check_file(args.pred),'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('Species\t') or not line:continue
			tmp = line.split('\t')
			if int(tmp[-2]) == 1:
				spe_list.append(tmp[0])
	gcf_spe_dic = {} # {gcf:spe}
	with gzip.open(check_file(args.classfy), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			if tmp[7] in spe_list:
				gcf_spe_dic.setdefault(tmp[0], tmp[7])
	tag_gcf_dic = {} # {tag:[gcf]}
	with gzip.open(check_file(args.database), 'rt') as IN:
		for line in IN:
			if line.startswith('>'):
				id = line.lstrip('>').split('|')[0]
			else:
				if id in gcf_spe_dic:
					tag_gcf_dic.setdefault(line.strip(), []).append(id)
	with gzip.open(args.output, 'wt') as OUT:
		for tag, gcf_list in tag_gcf_dic.items():
			tmp_spe_list = [] # {spe:None}
			for gcf in gcf_list:
				tmp_spe_list.append(gcf_spe_dic[gcf])
			if len(set(tmp_spe_list)) == 1:
				for gcf in set(gcf_list):
					OUT.write(">{}|1\n{}\n".format(gcf, tag))

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
