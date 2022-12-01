#!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, glob, joblib
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
import pandas as pd
import numpy as np
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/08/01 13:43:50'
__version__ = '1.0.0'
########################################### update ###############################################
"""
20220801 v1
单独预测功能
"""
############################################ main ##################################################
def report(level,info):
	date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	if level == "ERROR":
		sys.stderr.write("{0} - {1} - ERROR - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	elif level == "INFO":
		sys.stdout.write("{0} - {1} - INFO - {2}\n".format(date_now,os.path.basename(__file__),info))
	return()

def check_file(file):
	if os.path.exists(file):
		return(os.path.abspath(file))
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def pre_tra_data(df):
	sub_df = df.copy()
	sub_df.insert(5, 'Coverage_log_norm', np.log(sub_df['Sequenced_Tag_Num']/sub_df['Theoretical_Tag_Num']))
	sub_df.insert(6, 'G_Score_log_norm', np.log(sub_df['G_Score']))
	sub_df.insert(7, 'Tag_Depth_log_norm', np.log(sub_df['Sequenced_Reads_Num']/sub_df['Sequenced_Tag_Num']))
	sub_df.insert(8, 'Theoretical_Reads_DTR_log_norm', np.log((sub_df['Sequenced_Reads_Num']/sub_df['Sequenced_Tag_Num']*sub_df['Theoretical_Tag_Num'])/sub_df['Sequenced_Reads_Num'].sum()))
	return sub_df

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='test set (input file), qual or (qual + tag), allow to use *',dest='input',type=str,required=True)
	parser.add_argument('-m',help='module file, default RF_none_0238.v2.pkl',dest='module',type=str,default='{}/../tools/RF_none_0238.v2.pkl'.format(bindir))
	parser.add_argument('-y',help='output type: 1 for sub_spe, 2 for sub_spe to spe, 3 for sub_spe to spe to sub_spe, default 1',dest='type',choices=[1, 2, 3], type=int, default=1)
	parser.add_argument('-o',help='predicted result (output file)',dest='output',type=str,required=False)
	parser.add_argument('-t',help='the threshold for filtering through Gscore, default 5',dest='threshold',type=int,default=5)
	parser.add_argument('-s',help='species that need to be screened, use commas to join index, default human',dest='screened',type=str,default='human')

	args=parser.parse_args()
	feature_list = ['Coverage_log_norm', 'G_Score_log_norm', 'Tag_Depth_log_norm', 'Theoretical_Reads_DTR_log_norm']

	clf = joblib.load(args.module)

	file_list = glob.glob(args.input)
	if not file_list:
		report('ERROR', 'Error to find: {}'.format(args.input))

	capital_list = [chr(i) for i in range(65, 91)]
	for file in file_list:
		file = check_file(file)
		name = file.split('/')[-1]
		df = pd.read_csv(file, sep = '\t')
		if len(list(df)) == 1:
			df = pd.read_csv(file, skiprows=[0], sep = '\t')
		filter_s = 'Species != "' + '" & Species != "'.join(args.screened.split(',')) + '"'
		df = df.query(filter_s)
		df = df[(df['G_Score'] > args.threshold)].reset_index(drop=True)
		sub_df = df.iloc[:,[6, 8, 10, 7, 14]] # Species  Sequenced_Tag_Num  Sequenced_Reads_Num  Theoretical_Tag_Num  G_Score
		tra_df = pre_tra_data(sub_df)
		if args.type == 1:
			pred_list = clf.predict(tra_df[feature_list])
			prob_list = clf.predict_proba(tra_df[feature_list])[:,1]

		else:
			merge_spe_list = []
			merge_spe_dic = {}
			for i in df['Species']:
				if i[-1] in capital_list:
					merge_spe = '_'.join(i.split('_')[:-1])
					merge_spe_list.append(merge_spe)
					merge_spe_dic.setdefault(merge_spe, []).append(i)
				else:
					merge_spe_list.append(i)
					merge_spe_dic.setdefault(i, []).append(i)
			c_sub_df_tmp = df.copy()
			c_sub_df_tmp['Species'] = merge_spe_list
			c_sub_df = c_sub_df_tmp.groupby('Species').agg({'Sequenced_Tag_Num':'sum', 'Sequenced_Reads_Num':'sum', 'Theoretical_Tag_Num':'sum'}).reset_index(drop=False)
			c_sub_df['G_Score'] = np.sqrt(c_sub_df['Sequenced_Tag_Num'] * c_sub_df['Sequenced_Reads_Num'])
			c_tra_df = pre_tra_data(c_sub_df)
			spe_list = list(c_tra_df['Species'])
			pred_list = clf.predict(c_tra_df[feature_list])
			prob_list = clf.predict_proba(c_tra_df[feature_list])[:,1]
			if args.type == 2:
				tra_df = c_tra_df
				pred_list = pred_list
				prob_list = prob_list
			elif args.type == 3:
				pred_dic = dict(zip(spe_list, pred_list))
				prob_dic = dict(zip(spe_list, prob_list))
				pred_list, prob_list = [], []
				for s in merge_spe_list:
					pred_list.append(pred_dic[s])
					prob_list.append(prob_dic[s])

		if 'tag' in df:
			if args.type == 2:
				tag_sub_df = c_sub_df_tmp.groupby('Species').agg({'tag':'sum'}).reset_index(drop=False)
				tag_sub_df['tag'][tag_sub_df.tag > 1] = 1
				tra_df['Tag'] = tag_sub_df['tag']
			else:
				tra_df['Tag'] = df['tag']
		tra_df['Prediction'] = pred_list
		tra_df['Probability'] = prob_list
		if args.output:
			tra_df.to_csv(args.output, sep = '\t', index=None)
		else:
			tra_df.to_csv('pred_{}'.format(name), sep = '\t', index=None)

if __name__=="__main__":
	main()
