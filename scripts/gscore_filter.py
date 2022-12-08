#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3
########################################## import ################################################
import argparse, os, sys, re, random, glob
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/12/08 08:22:58'
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
	parser.add_argument('-g',help='gscore threshold',dest='gscore',type=int,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)
	I = check_file(args.input)
	n = 0
	with open(I,'r') as IN, open(args.output, 'w') as OUT:
		OUT.write('Species\tGscore\tPrediction\tProbability\n')
		for line in IN:
			line = line.strip()
			if line.startswith('#'):continue
			tmp = line.split('\t')
			if float(tmp[-1]) >= args.gscore:
				pred = 1
				n += 1
			else:
				pred = 0
			OUT.write('{spe}\t{gscore}\t{pred}\t{pred}\n'.format(spe = tmp[6], gscore = tmp[-1], pred = pred))
	if n == 0:
		report('ERROR', 'No microorganisms were found in this sample')

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)