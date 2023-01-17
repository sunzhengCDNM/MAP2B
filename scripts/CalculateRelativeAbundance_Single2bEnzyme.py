#!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip, marisa_trie
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from math import *
bindir = os.path.abspath(os.path.dirname(__file__))
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/11/20 02:31:28'
__version__ = '2.1.0'
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

def qual(outdir, smp, fasta):
	report('INFO', f"for sample {smp}")
	transtable = str.maketrans('ATGCN', 'TACGN')
	gcf_count_dic = {} # {gcf:{tag:seq_num}}
	tag_count_dic = {} # {tag:seq_num}
	# 读取fa，对tag进行去重复及统计
	with gzip.open(check_file(fasta), 'rt') as IN:
		for line in IN:
			if line.startswith('>'):continue
			seq_1 = line.strip()
			seq_2 = seq_1[::-1].translate(transtable)
			tag_count_dic[seq_1] = tag_count_dic.get(seq_1, 0) + 1
			tag_count_dic[seq_2] = tag_count_dic.get(seq_2, 0) + 1

	for tag, tag_num in tag_count_dic.items():
		try:
			for j in tag_gcf_trie[tag]:
				gcf = ''.join([str(i, 'utf-8') for i in j])
				tag_count_tmp_dic = gcf_count_dic.setdefault(gcf, {})
				tag_count_tmp_dic[tag] = tag_num
		except:
			continue

	tax_tag_count_dic = {} # {tax:tag_num}
	tax_reads_count_dic = {} # {tax:reads_num}
	with open(f"{outdir}/{smp}.{enzyme}.GCF_detected.xls", 'w') as OUT:
		for gcf, tag_count_dic in gcf_count_dic.items():
			tax = gcf_tax_dic[gcf]
			theo_tag_num = int(gcf_theo_tag_num_dic[gcf])
			dete_tag_num = len(tag_count_dic)
			tax_tag_count_dic.setdefault(tax, set()).update(set(tag_count_dic.keys()))
			for tag, tag_count in tag_count_dic.items():
				tmp_tax_reads_count_dic = tax_reads_count_dic.setdefault(tax, {})
				tmp_tax_reads_count_dic.setdefault(tag, tag_count)
			OUT.write('{}\t{}\t{}\t{}\t{}\n'.format(tax, gcf, theo_tag_num, dete_tag_num, round(dete_tag_num/theo_tag_num, 4)))
	with open(f"{outdir}/{smp}.{enzyme}.xls", 'w') as OUT:
		OUT.write('#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\tSequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\tSequenced_Reads_Num/Sequenced_Tag_Num\tSequenced_Tag_Num(depth>1)\tG_Score\n')
		for tax in tax_tag_count_dic:
			Theoretical_Tag_Num = float(tax_theo_tag_num_dic[tax])
			Sequenced_Tag_Num = len(tax_tag_count_dic[tax])
			Sequenced_Reads_Num = sum(tax_reads_count_dic[tax].values())
			OUT.write('{tax}\t{Theoretical_Tag_Num}\t{Sequenced_Tag_Num}\t{Percent}\t{Sequenced_Reads_Num}\t{Sequenced_Reads_Num__Theoretical_Tag_Num}\t{Sequenced_Reads_Num__Sequenced_Tag_Num}\tNA\t{G_Score}\n'.format(tax = tax, Theoretical_Tag_Num = Theoretical_Tag_Num, Sequenced_Tag_Num = Sequenced_Tag_Num, Percent = round(Sequenced_Tag_Num/Theoretical_Tag_Num, 4), Sequenced_Reads_Num = Sequenced_Reads_Num, Sequenced_Reads_Num__Theoretical_Tag_Num = round(Sequenced_Reads_Num/Theoretical_Tag_Num, 4), Sequenced_Reads_Num__Sequenced_Tag_Num = round(Sequenced_Reads_Num/Sequenced_Tag_Num, 4), G_Score = round(sqrt(Sequenced_Tag_Num*Sequenced_Reads_Num), 4)))

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-d',help='database prefix, like /path/CjePI.species.uniq',dest='database',type=str,required=True)
	parser.add_argument('-c',help='classify file, like abfh_classify_with_speciename.txt.gz',dest='classify',type=str,required=True)
	parser.add_argument('-l',help='sample list, sample <tab> enzyme.fa.gz',dest='list',type=str,required=True)
	parser.add_argument('-e',help='enzyme',dest='enzyme',type=str,required=True)
	parser.add_argument('-o',help='output dir',dest='output',type=str,required=True)
	parser.add_argument('-p',help='number of processes used',dest='processes',type=int,default=1)
	args=parser.parse_args()
	global enzyme, gcf_theo_tag_num_dic, tax_theo_tag_num_dic, gcf_tax_dic, tag_gcf_trie

# reading classify file
	report('INFO', 'Start reading classify file')
	gcf_tax_dic = {} # {gcf:tax}
	with gzip.open(check_file(args.classify), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			gcf_tax_dic.setdefault(tmp[0], '\t'.join(tmp[1:8]))
	report('INFO', 'End of reading classify file')

# reading database
	report('INFO', 'Start reading the database')
	db_file = check_file(args.database + '.marisa')
	stat_file = check_file(args.database + '.stat.xls')
	enzyme = args.enzyme
	tag_gcf_trie = marisa_trie.RecordTrie('8c').mmap(db_file)

# reading stat file
	gcf_theo_tag_num_dic = {} # {gcf: uniq_tag_num}
	tax_theo_tag_num_dic = {} # {tax: uniq_tag_num}
	with open(stat_file, 'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			gcf_theo_tag_num_dic.setdefault(tmp[0], tmp[8])
			tax_theo_tag_num_dic.setdefault('\t'.join(tmp[1:8]), tmp[9])
	report('INFO', 'End of reading database')

# qual
	report('INFO', 'Start qualitative analysis using {} processes'.format(args.processes))
	executor = ProcessPoolExecutor(args.processes)
	pool = []
	with open(check_file(args.list), 'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			outdir = check_dir('{}/{}'.format(args.output, tmp[0]))
			pool.append(executor.submit(qual, outdir, tmp[0], tmp[1]))
	executor.shutdown()
	for res in pool:
		res.result()
	report('INFO', 'End of qualitative analysis')

if __name__=="__main__":
	main()
