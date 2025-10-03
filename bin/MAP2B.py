#!python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
src_dir = os.path.abspath(os.path.dirname(__file__) + '/../scripts')
def_db_dir = os.path.abspath(os.path.dirname(__file__) + '/../database')
config_dir = os.path.abspath(os.path.dirname(__file__) + '/../config')

enzyme_dic = {1:'CspCI', 2:'AloI', 3:'BsaXI', 4:'BaeI',
      5:'BcgI', 6:'CjeI', 7:'PpiI', 8:'PsrI', 9:'BplI',
      10:'FalI', 11:'Bsp24I', 12:'HaeIV', 13:'CjePI',
      14:'Hin4I',15:'AlfI',16:'BslFI'}
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Zheng Sun, Liu Jiang'
__mail__ = 'spzsu@channing.harvard.edu, jiang.liu@oebiotech.com'
__date__ = '2025/10/02 20:03:47'
__version__ = '1.6'
############################################ main ##################################################
def report(level, info):
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

def check_dir(dir):
	dir = os.path.abspath(dir)
	if not os.path.exists(dir):
		os.system("mkdir -p {0}".format(dir))
		info = "mkdir: {0}".format(dir)
		report("INFO",info)
	return dir

def check_db(db, enzyme, source):
	db_lst = '{config_dir}/{source}.{enzyme}.database.list'.format(config_dir = config_dir, enzyme = enzyme, source = source)
	t = 0
	if not os.path.exists(db_lst):
		report('ERROR', 'Unable to find {db_lst}, please go to github (https://github.com/sunzhengCDNM/MAP2B/tree/master/config) to download'.format(db_lst = db_lst))
	else:
		with open(db_lst, 'r') as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('#') or not line:continue
				tmp = line.split('\t')
				db_file = '{}/{}'.format(db, tmp[0])
				if not os.path.exists(db_file):
					t = 1
	if t == 1:
		report('INFO', 'Downloading database, due to network reasons, may take a long time, please wait')
		exe_shell('python3 {src_dir}/DownloadDB.py -l {config_dir}/{source}.{enzyme}.database.list -d {db_dir}'.format(src_dir = src_dir, db_dir = db, config_dir = config_dir, enzyme = enzyme, source = source), 'DownloadDB')
	return

def exe_shell(cmd, desc=None):
	if desc:
		report('INFO', 'Running: {}'.format(desc))
	else:
		report('INFO', 'Running: {}'.format(cmd))
	code = os.system(cmd)
	if code == 0:
		return
	else:
		if desc:
			report('ERROR', 'failed to run: {}'.format(desc))
		else:
			report('ERROR', 'failed to run: {}'.format(cmd))

def mkdb(db_dir, enzyme, smp, quan_db, O, cutoff):
	exe_shell('python3 {src_dir}/creatQuanDB.py -d {db}/{enzyme}.species -c {db}/abfh_classify_with_speciename.txt.gz -p {O}/1.qual/{smp}/pred.result -o {quan_db}/{smp}.{enzyme} -e {enzyme} -s {cutoff} -n m'.format(src_dir = src_dir, db = db_dir, enzyme = enzyme, smp = smp, quan_db = quan_db, O = O, cutoff = cutoff), 'creatQuanDB')
	return

def check_data(data_file):
	data_dic = {}
	with open(data_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			data_dic.setdefault(tmp[0], [])
			for f in [check_file(i) for i in tmp[1:] if i]:
				data_dic.setdefault(tmp[0], []).append(f)
	for smp, data_list in data_dic.items():
			if len(data_list) > 2:
				report('ERROR', 'Sample {0} has more than 2 Reads, or sample {0} is named duplicate: {1}'.format(smp, ', '.join(data_list)))
			elif len(data_list) == 0:
				report('ERROR', 'Sample {} has no Reads'.format(smp))
	return data_dic

def cc_abd(db, abfh, enzyme_smp_file, O, processes, enzyme):
	exe_shell('python3 {src_dir}/CalculateRelativeAbundance_Single2bEnzyme.py -d {db} -c {abfh} -l {enzyme_smp_file} -o {O} -p {processes} -e {enzyme}'.format(src_dir = src_dir, db = db, abfh = abfh, enzyme_smp_file = enzyme_smp_file, O = O, processes = processes, enzyme = enzyme), 'CalculateRelativeAbundance_Single2bEnzyme')

def extra_tag(reads, enzyme, enzyme_dir, smp):
	if len(reads) == 2:
		exe_shell('perl {src_dir}/2bRADExtraction.pl -i {reads} -t 2 -s {enzyme} -od {enzyme_dir} -op {smp}_1 -qc no'.format(src_dir = src_dir, reads = reads[0], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
		exe_shell('perl {src_dir}/2bRADExtraction.pl -i {reads} -t 2 -s {enzyme} -od {enzyme_dir} -op {smp}_2 -qc no'.format(src_dir = src_dir, reads = reads[1], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
		exe_shell('cat {enzyme_dir}/{smp}_1.{enzyme}.fa.gz {enzyme_dir}/{smp}_2.{enzyme}.fa.gz >{enzyme_dir}/{smp}.{enzyme}.fa.gz && rm {enzyme_dir}/{smp}_1.{enzyme}.fa.gz {enzyme_dir}/{smp}_2.{enzyme}.fa.gz'.format(enzyme = enzyme_dic[enzyme], enzyme_dir = enzyme_dir, smp = smp), 'mergePEReads')
	else:
		exe_shell('perl {src_dir}/2bRADExtraction.pl -i {reads} -t 2 -s {enzyme} -od {enzyme_dir} -op {smp} -qc no'.format(src_dir = src_dir, reads = reads[0], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
	return

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\nlast update:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. like \n \
	sample <tab> shotgun.1.fq(.gz) (<tab> shotgun.2.fq.gz)',dest='input',type=str,required=True)
	parser.add_argument('-o',help='Output directory, default {}/MAP2B_result'.format(os.getcwd()),dest='output',type=str,default='{}/MAP2B_result'.format(os.getcwd()))
	parser.add_argument('-e',help='Enzyme, choose from 5(BcgI) or 13(CjePI), default 13',dest='enzyme',type=int,choices=[5, 13],default=13)
#	parser.add_argument('-e',help='enzyme, default 13 for CjePI, choose from\n \
#	[1]CspCI  [5]BcgI  [9]BplI     [13]CjePI  [17]AllEnzyme\n \
#	[2]AloI   [6]CjeI  [10]FalI    [14]Hin4I\n \
#	[3]BsaXI  [7]PpiI  [11]Bsp24I  [15]AlfI\n \
#	[4]BaeI   [8]PsrI  [12]HaeIV   [16]BslFI',dest='enzyme',type=int,default=13)
	parser.add_argument('-s',help='Data source, choose from GTDB or RefSeq, default GTDB',dest='source',type=str,choices=['GTDB', 'RefSeq'],default='GTDB')
	parser.add_argument('-d',help='Database path for MAP2B pipeline, default {}'.format(def_db_dir),dest='database',type=str,default=def_db_dir)
	parser.add_argument('-p',help='Number of processes, note that more threads may require more memory, default 1',dest='processes',type=int,default=1)
	parser.add_argument('-g',help='Using G score as the threshold for species identification, -g 5 is recommended. Enabling G score will automatically shutdown false positive recognition model, default none',dest='gscore',type=int,required=False)
	parser.add_argument('-c',help='cut off for database, default 30000',dest='cutoff',type=int,default=30000)

	args=parser.parse_args()
	enzyme_id = args.enzyme
#	enzyme_id = 13
	enzyme = enzyme_dic[enzyme_id]
	db_dir = check_dir((args.database + '/' + args.source))
	# check database
	check_db(db_dir, enzyme, args.source)
	O = check_dir(args.output)
	# prepare data
	data_dic = check_data(check_file(args.input))
	quan_smp_set = set(data_dic.keys())
	enzyme_dir = check_dir(O + '/0.dige')

	done_file = O + '/0.dige/done'
	enzyme_smp_file = O + '/0.dige/enzyme_smp.list'
	if os.path.exists(done_file):
		report('INFO', 'The data digestion has been completed, go to the next step')
	else:
		report('INFO', 'Start of data digestion')
		check_dir(O + '/0.dige')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		with open(enzyme_smp_file, 'w') as OUT:
			for smp, reads in data_dic.items():
				OUT.write('{smp}\t{enzyme_dir}/{smp}/{smp}.{enzyme}.fa.gz\n'.format(enzyme_dir = enzyme_dir, smp = smp, enzyme = enzyme))
				pool.append(executor.submit(extra_tag, reads, enzyme_id, (enzyme_dir + '/' + smp), smp))
		executor.shutdown()
		for res in pool:
			res.result()
		exe_shell('touch {}'.format(done_file), 'dige_done')

	done_file = O + '/1.qual/done_q'
	if os.path.exists(done_file):
		report('INFO', 'The qualitative analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start qualitative analysis')
		check_dir(O + '/1.qual')
		cc_abd('{}/{}.species.uniq'.format(db_dir, enzyme), '{}/abfh_classify_with_speciename.txt.gz'.format(db_dir), enzyme_smp_file, O + '/1.qual', args.processes, enzyme)
		exe_shell('touch {}'.format(done_file), '')

	done_file = O + '/1.qual/done_p'
	none_spe_smp_list = []
	if os.path.exists(done_file):
		report('INFO', 'The predictive analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start predictive analysis')
		check_dir(O + '/1.qual')
		for smp in data_dic.keys():
			try:
				if args.gscore:
					exe_shell('python3 {src_dir}/gscore_filter.py -i {O}/1.qual/{smp}/{smp}.{enzyme}.xls -o {O}/1.qual/{smp}/pred.result -g {gscore}'.format(src_dir = src_dir, O = O, smp = smp, enzyme= enzyme, gscore = args.gscore), 'gscore_filter')
				else:
					exe_shell('python3 {src_dir}/MAP2B_ML.py -i {O}/1.qual/{smp}/{smp}.{enzyme}.xls -m {config_dir}/RF_none_0238.v2.pkl -o {O}/1.qual/{smp}/pred.result'.format(src_dir = src_dir, O = O, smp = smp, enzyme= enzyme, config_dir = config_dir), 'MAP2B_ML')
			except:
				report('INFO', 'If an error is reported there, it may be due to the absence of microorganisms, so please ignore it')
				none_spe_smp_list.append(smp)
			tmp_p_list = []
			with open('{O}/1.qual/{smp}/pred.result'.format(O = O, smp = smp), 'r') as IN:
				for line in IN:
					tmp = line.split()
					if tmp[0] == 'Species':continue
					if int(tmp[-2]) == 1:
						tmp_p_list.append(1)
			if not tmp_p_list:
				none_spe_smp_list.append(smp)
		quan_smp_set = set(data_dic.keys()) - set(none_spe_smp_list)
		if not quan_smp_set:
			report('ERROR', 'No microorganisms were found in all samples')
		exe_shell('touch {}'.format(done_file), 'filter_done')

	done_file = O + '/2.mkdb/done'
	if os.path.exists(done_file):
		report('INFO', 'The quantitative database building has been completed, go to the next step')
	else:
		report('INFO', 'Start quantitative database building')
		check_dir(O + '/2.mkdb')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		for smp in quan_smp_set:
			quan_db = '{O}/2.mkdb/{smp}'.format(O = O, smp = smp)
			check_dir(quan_db)
			with open('{}/reads.list'.format(quan_db), 'w') as OUT:
				OUT.write('{smp}\t{enzyme_dir}/{smp}/{smp}.{enzyme}.fa.gz\n'.format(enzyme_dir = enzyme_dir, smp = smp, enzyme = enzyme))
			pool.append(executor.submit(mkdb, db_dir, enzyme, smp, quan_db, O, args.cutoff))
		executor.shutdown()
		for res in pool:
			res.result()
		exe_shell('touch {}'.format(done_file), 'mkdb_done')

	done_file = O + '/3.quan/done'
	abd_list = '{}/3.quan/abd.list'.format(O)
	if os.path.exists(done_file):
		report('INFO', 'The quantitative analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start quantitative analysis')
		check_dir(O + '/3.quan')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		for smp in quan_smp_set:
			pool.append(executor.submit(cc_abd, '{}/2.mkdb/{}/{}.{}'.format(O, smp, smp, enzyme), '{}/abfh_classify_with_speciename.txt.gz'.format(db_dir), '{}/2.mkdb/{}/reads.list'.format(O, smp), O + '/3.quan', 1, enzyme))
		executor.shutdown()
		for res in pool:
			res.result()
		with open(abd_list, 'w') as OUT:
			for smp in quan_smp_set: 
				OUT.write('{}\t{}/3.quan/{}/{}.{}.xls\n'.format(smp, O, smp, smp, enzyme))
		exe_shell('touch {}'.format(done_file), 'quan_done')

	done_file = O + '/all_done'
	if os.path.exists(done_file):
		pass
	else:
		exe_shell('perl {}/MergeProfilesFromMultipleSamples.pl -l {} -o {} -p Abundance -m -c'.format(src_dir, abd_list, O), 'MergeProfilesFromMultipleSamples')
		exe_shell('python3 {}/MergeCoverageFromMultipleSamples.py -i {} -o {}/Coverage.xls'.format(src_dir, abd_list, O), 'MergeCoverageFromMultipleSamples')
		exe_shell('touch {}'.format(done_file), 'all_done')
	report('INFO', 'Congratulations, all work has been completed')

if __name__=="__main__":
	main()
