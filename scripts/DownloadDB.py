#!python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, hashlib
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2023/01/16 22:12:02'
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

def get_file_md5(file_path):
	with open(file_path, 'rb') as f:
		md5obj = hashlib.md5()
		md5obj.update(f.read())
		_hash = md5obj.hexdigest()
		return str(_hash)

def check_dbfile(local_file, md5):
	if not os.path.exists(local_file):
		return False
	else:
		if get_file_md5(local_file) == md5:
			report('INFO', 'File {} checksum passed'.format(local_file))
			return True
		else:
			return False

def downloaded_dbfile(url, local_file, md5):
	for i in range(3):
		try:
			os.system('wget --user-agent="Mozilla/5.0" -t 3 -O {local_file} {url}'.format(local_file = local_file, url = url))
			if check_dbfile(local_file, md5):
				return True
			else:
				continue
		except:
			pass
	return False

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-l',help='Database file list',dest='dblist',type=str,required=True)
	parser.add_argument('-d',help='Database path',dest='dbpath',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)
	I = check_file(args.dblist)
	O = check_dir(args.dbpath)
	with open(I,'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			name, url, md5 = line.split('\t')
			local_file = '{}/{}'.format(O, name)
			if check_dbfile(local_file, md5):
				continue
			else:
				if downloaded_dbfile(url, local_file, md5):
					continue
				else:
					report('ERROR', '{} download failed, please check the network'.format(url))
	report('INFO', 'All database files have been downloaded, please start the analysis')

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
