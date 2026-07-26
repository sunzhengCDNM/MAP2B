#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3

"""
数据库结构为
key： 8 位 genome_ID + 4 位 tag copy num(最大到 9999)
value: tag

v2.0 20230220
1. 修复基因组电子酶切时跳过重叠 tag 的 bug。

v2.1 20230223
1. 使用正则表达式的“正向前瞻”和“零宽度断言”功能进行重叠 tag 的识别。

"""
########################################## import ################################################
import argparse, os, sys, marisa_trie, gzip, re, collections, math
from datetime import datetime
############################################ ___ #################################################
__doc__ = 'This program is used for MAP2B database building'
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2023/01/08 23:33:25'
__version__ = '2.3'
############################################ main ##################################################
enzyme_pattern_dic = {
					'AlfI':r'(?=([AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}))',
					'AloI':r'(?=([AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}))',
					'BaeI':r'(?=([AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}))',
					'BcgI':r'(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))',
					'BplI':r'(?=([AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}))',
					'BsaXI':r'(?=([AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}))',
					'BslFI':r'(?=([AGCT]{6}GGGAC[AGCT]{14}))',
					'Bsp24I':r'(?=([AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}))',
					'CjeI':r'(?=([AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}))',
					'CjePI':r'(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))',
					'CspCI':r'(?=([AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}))',
					'FalI':r'(?=([AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}))',
					'HaeIV':r'(?=([AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}))',
					'Hin4I':r'(?=([AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}))',
					'PpiI':r'(?=([AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}))',
					'PsrI':r'(?=([AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}))',
					}

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
		return file
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def save_db(out, n, Tag_lst, ID_lst):
	trie = marisa_trie.BytesTrie(zip(ID_lst, Tag_lst))
	trie.save('{}.marisa{}'.format(out, n))

def read_fa(genome):
	seq_list = []
	seq = ''
	try:
		with gzip.open(genome, 'rt') as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('>'):
					seq_list.append(seq)
					seq = ''
				else:
					seq += line.upper()
		seq_list.append(seq)
	except EOFError:
		report("INFO", "警告：基因组文件损坏或不完整，跳过: {0}".format(genome))
		return []
	except OSError as e:
		report("INFO", "警告：读取基因组文件失败，跳过 {0}: {1}".format(genome, str(e)))
		return []
	return [seq for seq in seq_list if seq]

def extraction(genome, enzyme_pattern_list):
	tag_num_dic = collections.defaultdict(int)
	for seq in read_fa(genome):
		for enzyme_pattern in enzyme_pattern_list:
			for tag in re.findall(enzyme_pattern, seq):
				tag_num_dic[tag.zfill(40)] += 1
#		for enzyme_pattern in enzyme_pattern_list:
			for tag in re.findall(enzyme_pattern, seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))):
				tag_num_dic[tag.zfill(40)] += 1
	return tag_num_dic

def progress_bar(p):
	sys.stdout.write('\r{}{} {}%'.format('#' * p, '-' * (50 - p), p * 2))
	sys.stdout.flush()

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='abfh_classify_with_speciename.txt.gz 文件，第一列为 00000000 开始的连续、顺序编号',dest='classfy',type=str,required=True)
	parser.add_argument('-e',help='酶切位点组合，可选 AlfI/AloI/BaeI/BcgI/BplI/BsaXI/BslFI/Bsp24I/CjeI/CjePI/CspCI/FalI/HaeIV/Hin4I/PpiI/PsrI 中的一个或多个，多个酶之间用逗号间隔，如果选择全部酶可直接提供 all',dest='enzyme',type=str, required=True)
	parser.add_argument('-s',help='数据库拆分大小，表示每个子数据库中包含的基因组数，建议 30000，需要注意的是该值越大，建库所需内存越大',dest='size',type=int,required=True)
	parser.add_argument('-o',help='输出数据库前缀, 例如 BcgI.species',dest='output',type=str,required=True)
#	parser.add_argument('-f',help='需要过滤掉的低频 tag（最少 N 个基因组中唯一出现的 tag），建议 100, 无需过滤时省略该参数或者设置为 0 即可',dest='frequency',type=int,required=False)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)

	n, m, k = 0, 0, 0
	if args.enzyme == 'all':
		enzyme_pattern_list = list(enzyme_pattern_dic.values())
	else:
		enzyme_pattern_list = [enzyme_pattern_dic[enzyme] for enzyme in args.enzyme.split(',')]
	thd = int(args.size)

	# 从 abfh_classify_with_speciename.txt.gz 构建全数据库
	# 读取 abfh 文件
	class_dic = {} # {spe:{ID:genome}}
	with gzip.open(check_file(args.classfy), 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			tmp_dic = class_dic.setdefault(tmp[7], {})
			tmp_dic.setdefault(tmp[0], tmp[-1])
			m += 1

	# 统计 tag copy 数
	x = thd
	marisa_ID_lst, marisa_Tag_lst = [], []
	for spe, genome_dic in class_dic.items():
		sID_tag_dic = collections.defaultdict(list)  # {ID_tag_num:[tag]}
		for ID, genome in genome_dic.items():
			tag_num_dic = extraction(genome, enzyme_pattern_list)
			if not tag_num_dic:
				# 进度条（即使无数据也更新）
				n += 1
				p = int(n / m * 50)
				if p != k:
					progress_bar(p)
					k = p
				continue
			for tag, tag_num in tag_num_dic.items():
				sID = ID + str(tag_num).zfill(4)
				sID_tag_dic[sID].append(tag)
			# 进度条
			n += 1
			p = int(n / m * 50)
			if p == k:continue
			progress_bar(p)
			k = p

#		# 过滤低频 tag
#		tag_ID_dic = collections.defaultdict(list)  # {tag:[ID]}
#		ID_tag_dic = collections.defaultdict(list)  # {ID:[tag]}
#		uniq_tag_lst = []
#		if args.frequency:
#			if len(genome_dic) > args.frequency:
#				for tag, ID_lst in tag_ID_dic.items():
#					if len(ID_lst) == 1:
#						uniq_tag_lst.append(tag)
#		for tag in uniq_tag_lst:
#			tag_ID_dic.pop(tag)
#		for tag, ID_lst in tag_ID_dic.items(): 
#			for ID in ID_lst:
#				ID_tag_dic[ID].append(tag)

		# 分段存储数据库
		# v2.2 修改：支持非零起始编号
		# 仅在全局列表为空时（即第一个有数据的物种），调整x到正确的起始阈值
		if not marisa_Tag_lst and sID_tag_dic:
			first_ID = int(sorted(sID_tag_dic.keys())[0][:8])
			if first_ID >= x:
				x = ((first_ID - 1) // thd + 1) * thd
		for sID in sorted(sID_tag_dic.keys()):
			ID = int(sID[:8])
			if ID >= x:
				if marisa_Tag_lst:
					save_db(args.output, x, marisa_Tag_lst, marisa_ID_lst)
					x += thd
					marisa_ID_lst, marisa_Tag_lst = [], []
				else:
					# 空列表时只更新x，不保存空数据库
					x += thd
			tag_lst = [tag.encode('utf-8') for tag in sID_tag_dic[sID]]
			marisa_Tag_lst += tag_lst
			marisa_ID_lst += [sID] * len(tag_lst)
	if marisa_Tag_lst:
		save_db(args.output, x, marisa_Tag_lst, marisa_ID_lst)

if __name__=="__main__":
	main()
	sys.stdout.write('\n')
	info = "finish!"
	report("INFO",info)
