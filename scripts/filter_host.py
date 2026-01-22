#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3
########################################## import ################################################
import argparse, os, sys, re, random, glob
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append('/data/USER/liujiang/script/lib')
#import 
#import pandas as pd
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Zheng & Jiang'
__mail__ = 'spzsu@channing.harvard.edu'
__date__ = '2021/01/22 15:56:33'
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
	parser.add_argument('-i',help='input file',dest='input',type=str,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	parser.add_argument('-n',help='quantitative column, default 8',dest='n',type=int,default=8)
	parser.add_argument('-s',help='host name, human/mus, default human',dest='host',type=str,default='human')
	args=parser.parse_args()
	info = "runing..."
#	report("INFO",info)
	n = args.n - 1
	info_list = []
	I = check_file(args.input)
	none_index = []
	with open(I,'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#'):
				header = line
			else:
				tmp = line.split('\t')
				if tmp[0] == args.host:
					trans_list = []
					for i in tmp[n:]:
						try:
							trans_list.append((1/(1-float(i))))
						except:
							trans_list.append(0)
							none_index.append(tmp[n:].index(i))
				else:
					info_list.append(tmp)
#	print(header)
	with open(args.output,'w') as OUT:
#		h = '\t'.join(['_'.join(header[:n])] + header[n:]) + '\n'
#		OUT.write(h.strip('#'))
		OUT.write(header + "\n")
		for i in info_list:
			abun_list = [float(j) for j in i[n:]]
#			print(abun_list)
			try:
				new_abun_list = [a*b for a,b in zip(abun_list, trans_list)]
			except:
				new_abun_list = [a*b for a,b in zip(abun_list, [1]*len(abun_list))]
#			new_info_list = ['_'.join(i[:n])] + [str(k) for k in new_abun_list]
			new_info_list = i[:n] + [str(k) for k in new_abun_list]
			OUT.write('\t'.join(new_info_list) + '\n')
		if none_index:
			a = ["NONE"] * n
			b = ["0"] * (len(i) - n)
			for j in none_index:
				b[j] = "1"
			OUT.write('\t'.join(a + b) + '\n')

if __name__=="__main__":
	main()
	info = "finish!"
#	report("INFO",info)
