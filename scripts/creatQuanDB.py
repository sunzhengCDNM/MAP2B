 #!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip, marisa_trie, collections
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
	return

def check_file(file):
	if os.path.exists(file):
		return os.path.abspath(file)
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def ID_subdb(thd, ID_lst):
	subdb_ID_dic = {} # {classify:[ID]}
	for ID in sorted(ID_lst):
		subdb_ID_dic.setdefault((((int(ID) // thd) + 1) * thd), []).append(ID)
	return subdb_ID_dic

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-d',help='database prefix, like CjePI.species',dest='database',type=str,required=True)
	parser.add_argument('-c',help='classfy file',dest='classfy',type=str,required=True)
	parser.add_argument('-p',help='pred file',dest='pred',type=str,required=True)
	parser.add_argument('-s',help='the size of the sub-db when the master database is built',dest='size',type=int,required=True)
	parser.add_argument('-o',help='out_database prefix',dest='output',type=str,required=True)
	parser.add_argument('-e',help='enzyme, choose from BcgI and CjePI',dest='enzyme',choices=['BcgI', 'CjePI'],type=str,required=True)
	parser.add_argument('-n',help='copy number for tag, choos from s(ingle) and m(ultiple)',dest='copy',choices=['s', 'm'],type=str,required=True)
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
	ID_spe_dic = {} # {ID:spe}
	with gzip.open(check_file(args.classfy), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			if tmp[7] in spe_lst:
				ID_spe_dic.setdefault(tmp[0], '\t'.join(tmp[1:8]))
	tag_sID_dic = {} # {tag:[sID]}
	subdb_ID_dic = ID_subdb(args.size, ID_spe_dic.keys())
	for subdb_num, ID_lst in subdb_ID_dic.items():
		marisa = '{}.marisa{}'.format(args.database, subdb_num)
		trie = marisa_trie.RecordTrie(fmt).mmap(marisa)
		for ID in ID_lst:
			for sID in set(trie.keys(ID)):
				for tag in trie[sID]:
					tag_sID_dic.setdefault(''.join([str(i, 'utf-8') for i in tag]), []).append(sID)

	ID_lst_trie, uniq_tag_lst_trie= [], []
	ID_theo_tag_num_dic = collections.defaultdict(int) # {ID: uniq_tag_num}
	if args.copy == 'm':
		for tag, sID_lst in tag_sID_dic.items():
			tmp_spe_lst = [] # [spe]
			for sID in sID_lst:
				tmp_spe_lst.append(ID_spe_dic[sID[:8]])
			if len(set(tmp_spe_lst)) == 1:
				tmp_ID_lst_trie = []
				for sID in sID_lst:
					tmp_ID_lst_trie.append(sID[:8])
					ID_theo_tag_num_dic[sID[:8]] += int(sID[8:])
				ID_lst_trie += tmp_ID_lst_trie
				uniq_tag_lst_trie += [tag] * len(tmp_ID_lst_trie)
	elif args.copy == 's':
		for tag, sID_lst in tag_sID_dic.items():
			tmp_spe_lst = [] # [spe]
			for sID in sID_lst:
				tmp_spe_lst.append(ID_spe_dic[sID[:8]])
			if len(set(tmp_spe_lst)) == 1:
				tmp_ID_lst_trie = [sID[:8] for sID in sID_lst if sID[8:] == '0001']
				for ID in tmp_ID_lst_trie:
					ID_theo_tag_num_dic[ID] += 1
				ID_lst_trie += tmp_ID_lst_trie
				uniq_tag_lst_trie += [tag] * len(tmp_ID_lst_trie)

	del tag_sID_dic
	trie = marisa_trie.BytesTrie(zip(uniq_tag_lst_trie, [ID.encode('utf-8') for ID in ID_lst_trie]))
	del uniq_tag_lst_trie, ID_lst_trie
	trie.save(marisa_file)
	del trie

## 统计文件
	spe_theo_tag_num_dic = {} # {spe: uniq_tag_num}
	for ID, sum_t in ID_theo_tag_num_dic.items():
		spe_theo_tag_num_dic.setdefault(ID_spe_dic[ID], []).append(sum_t)
	for spe, t in spe_theo_tag_num_dic.items():
		spe_theo_tag_num_dic[spe] = round(sum(t)/len(t), 4)
	with open(stat_file, 'w') as OUT:
		for ID in sorted(ID_theo_tag_num_dic.keys()):
			OUT.write('{}\t{}\t{}\t{}\n'.format(ID, ID_spe_dic[ID], ID_theo_tag_num_dic[ID], spe_theo_tag_num_dic[ID_spe_dic[ID]]))

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
