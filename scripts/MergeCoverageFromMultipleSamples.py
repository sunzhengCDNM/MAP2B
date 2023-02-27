import argparse, os, sys, re, random, glob
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))

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
	parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-i',help='input file',dest='input',type=str,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
#	report("INFO",info)
	info_dic = {}
	I = check_file(args.input)
	with open(I,'r') as IN:
		for line1 in IN:
			smp, abd = line1.split()
			with open(abd, 'r') as ABD:
				for line2 in ABD:
					if line2.startswith('#'):continue
					tmp = line2.strip().split('\t')
					tmp_dic = info_dic.setdefault(smp, {})
					tmp_dic.setdefault('\t'.join(tmp[:7]), tmp[9])
	all_spe_list = []
	for smp, spe_cov_dic in info_dic.items():
		all_spe_list += spe_cov_dic.keys()
	with open(args.output, 'w') as OUT:
		smp_list = info_dic.keys()
		OUT.write('#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t{}\n'.format('\t'.join(smp_list)))
		for spe in sorted(list(set(all_spe_list))):
			cov_list = []
			for smp in smp_list:
				try:
					if float(info_dic[smp][spe]) > 1:
						cov = '1'
					else:
						cov = info_dic[smp][spe]
					cov_list.append(cov)
				except:
					cov_list.append('0')
			OUT.write('{}\t{}\n'.format(spe, '\t'.join(cov_list)))




		

if __name__=="__main__":
	main()
	info = "finish!"
#	report("INFO",info)
