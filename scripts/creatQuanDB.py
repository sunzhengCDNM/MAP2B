#!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip, marisa_trie
from datetime import datetime
############################################ ___ #################################################
__doc__ = '本程序用于MAP2B流程的定量建库'
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

def gcf_classify(thd, gcf_lst):
	gcf_classify_dic = {} # {classify:[gcf]}
	for gcf in sorted(gcf_lst):
		gcf_classify_dic.setdefault((((int(gcf) // thd) + 1) * thd), []).append(gcf)
	return gcf_classify_dic

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-d',help='database prefix, like CjePI.species',dest='database',type=str,required=True)
	parser.add_argument('-c',help='classfy file',dest='classfy',type=str,required=True)
	parser.add_argument('-p',help='pred file',dest='pred',type=str,required=True)
	parser.add_argument('-t',help='threshold for the database building',dest='threshold',type=int,required=True)
	parser.add_argument('-o',help='out_database prefix',dest='output',type=str,required=True)
	parser.add_argument('-e',help='enzyme, choose from BcgI and CjePI',dest='enzyme',type=str,default='CjePI')
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)
	spe_lst = []
	enzyme_dic = {'CjePI':27, 'BcgI':32}
	marisa_file = '{}.marisa'.format(args.output)
	stat_file = '{}.stat.xls'.format(args.output)
	fmt = '{}c'.format(enzyme_dic[args.enzyme])

	with open(check_file(args.pred),'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('Species\t') or not line:continue
			tmp = line.split('\t')
			if int(tmp[-2]) == 1:
				spe_lst.append(tmp[0])
	gcf_spe_dic = {} # {gcf:spe}
	with gzip.open(check_file(args.classfy), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			if tmp[7] in spe_lst:
				gcf_spe_dic.setdefault(tmp[0], '\t'.join(tmp[1:8]))
	tag_gcf_dic = {} # {tag:[gcf]}
	gcf_classify_dic = gcf_classify(args.threshold, gcf_spe_dic.keys())
	for classify, gcf_lst in gcf_classify_dic.items():
		marisa = '{}.marisa{}'.format(args.database, classify)
		trie = marisa_trie.RecordTrie(fmt).mmap(marisa)
		for gcf in gcf_lst:
			for tag in trie[gcf]: 
				tag_gcf_dic.setdefault(''.join([str(i, 'utf-8') for i in tag]), []).append(gcf)
	gcf_lst_trie, uniq_tag_lst_trie= [], []
	for tag, gcf_lst in tag_gcf_dic.items():
		tmp_spe_lst = [] # [spe]
		for gcf in gcf_lst:
			tmp_spe_lst.append(gcf_spe_dic[gcf])
		if len(set(tmp_spe_lst)) == 1:
			for gcf in set(gcf_lst):
				uniq_tag_lst_trie.append(tag)
				gcf_lst_trie.append(gcf)
	del tag_gcf_dic
	trie = marisa_trie.BytesTrie(zip(uniq_tag_lst_trie, [gcf.encode('utf-8') for gcf in gcf_lst_trie]))
	trie.save(marisa_file)
	trie = marisa_trie.BytesTrie(zip(gcf_lst_trie, [tag.encode('utf-8') for tag in uniq_tag_lst_trie]))
	del uniq_tag_lst_trie

## 统计文件
	gcf_theo_tag_num_dic = {} # {gcf: uniq_tag_num}
	spe_theo_tag_num_dic = {} # {spe: uniq_tag_num}
	for gcf in set(gcf_lst_trie):
		gcf_theo_tag_num_dic[gcf] = len(trie[gcf])
	del gcf_lst_trie, trie
	for gcf, sum_t in gcf_theo_tag_num_dic.items():
		spe_theo_tag_num_dic.setdefault(gcf_spe_dic[gcf], []).append(sum_t)
	for spe, t in spe_theo_tag_num_dic.items():
		spe_theo_tag_num_dic[spe] = round(sum(t)/len(t), 4)
	with open(stat_file, 'w') as OUT:
		for gcf, gcf_uniq_count in gcf_theo_tag_num_dic.items():
			OUT.write('{}\t{}\t{}\t{}\n'.format(gcf, gcf_spe_dic[gcf], gcf_theo_tag_num_dic[gcf], spe_theo_tag_num_dic[gcf_spe_dic[gcf]]))

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
