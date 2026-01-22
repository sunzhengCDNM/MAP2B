#!/data/software/install/miniconda3/envs/python.3.7.0/bin//home/liujiang/software/install/mamba/envs/map2b/bin/python3

"""
数据库结构为
key： 酶切出来的序列
value: 序列id__序列方向__位置索引（序列位置索引从 0 开始，反向序列的位置从最后一个字符往前计算）

"""
########################################## import ################################################
import argparse, os, sys, marisa_trie, gzip, re, collections, math
from datetime import datetime
############################################ ___ #################################################
__doc__ = '多线程、多存储结构的电子酶切脚本'
__author__ = 'Zheng and Jiang'
__mail__ = 'spzsu@channing.harvard.edu'
__date__ = '2024/06/01 23:33:25'
__version__ = '1.1'  # 版本更新
############################################ main ##################################################
enzyme_pattern_dic = {
    'BcgI': [32, r'(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))'],
    'CjePI': [27, r'(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))'],
    'BsaXI': [27, r'(?=([AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}))'],
    'CspCI': [32, r'(?=([AGCT]{10}CAA[AGCT]{5}GTGG[AGCT]{10}))']
}

def report(level, info):
    date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if level == "ERROR":
        sys.stderr.write(f"{date_now} - {os.path.basename(__file__)} - ERROR - {info}\n")
        sys.exit(1)
    elif level == "INFO":
        sys.stdout.write(f"{date_now} - {os.path.basename(__file__)} - INFO - {info}\n")
    elif level == "DEBUG":
        sys.stdout.write(f"{date_now} - {os.path.basename(__file__)} - DEBUG - {info}\n")
        sys.exit(1)
    return

def check_file(file):
    if os.path.exists(file):
        return file
    else:
        info = f"file does not exist: {file}"
        report("ERROR", info)

def read_cram(cram, fasta):
    import pysam
    with pysam.AlignmentFile(cram, mode="rc", reference_filename=fasta, require_index=True) as IN:
        for read in IN.fetch():
            yield read.query_name, read.seq

def open_file(sequence_file):
    if sequence_file.endswith('.gz'):
        return gzip.open(sequence_file, 'rt')
    else:
        return open(sequence_file, 'r')

def read_fafq(sequence_file):
    with open_file(sequence_file) as IN:
        try:
            first_line = IN.readline().strip()
        except UnicodeDecodeError:
            report('ERROR', '请给定 fa 或者 fq 格式的序列文件，不要提供 marisa 或者其他二进制的格式！')
    if first_line.startswith('>'):
        id = None
        seq = ''
        with open_file(sequence_file) as IN:
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    if id is not None:
                        yield id, seq
                    id = line.lstrip('>').split()[0]
                    seq = ''
                else:
                    seq += line.upper()
            if id is not None:
                yield id, seq
    elif first_line.startswith('@'):
        with open_file(sequence_file) as IN:
            while True:
                id_line = IN.readline()
                if not id_line:
                    break  # End of file
                id = id_line.split()[0]
                seq = IN.readline().strip()
                _ = IN.readline()  # Skip the separator line
                _ = IN.readline()  # Skip the quality score line
                yield id, seq
    else:
        report('ERROR', 'Unknown file format (not fastq or fasta)!')

def extraction(sequence_file, value_len, enzyme_pattern, genome=None):
    ori_seq_ct = 0
    dig_seq_ct = 0
    trans_table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    
    if genome and sequence_file.split('.')[-1] in ['cram', 'CRAM']:
        for id, seq in read_cram(sequence_file, genome):
            ori_seq_ct += 1
            # 正向序列匹配
            for matches in re.finditer(enzyme_pattern, seq):
                dig_seq_ct += 1
                new_id = f'{id}_{"+"}_{(matches.start() + 1)}'.rjust(value_len, '.')[-value_len:]
                yield matches.group(1), new_id
            # 反向互补序列匹配
            reversed_complement = seq[::-1].translate(trans_table).upper()
            for matches in re.finditer(enzyme_pattern, reversed_complement):
                dig_seq_ct += 1
                new_id = f'{id}_{"-"}_{(matches.start() + 1)}'.rjust(value_len, '.')[-value_len:]
                yield matches.group(1), new_id
    elif sequence_file.split('.')[-1] in ['cram', 'CRAM']:
        report('ERROR', 'give -g please')
    else:
        for id, seq in read_fafq(sequence_file):
            ori_seq_ct += 1
            # 正向序列匹配
            for matches in re.finditer(enzyme_pattern, seq):
                dig_seq_ct += 1
                new_id = f'{id}_{"+"}_{(matches.start() + 1)}'.rjust(value_len, '.')[-value_len:]
                yield matches.group(1), new_id
            # 反向互补序列匹配
            reversed_complement = seq[::-1].translate(trans_table).upper()
            for matches in re.finditer(enzyme_pattern, reversed_complement):
                dig_seq_ct += 1
                new_id = f'{id}_{"-"}_{(matches.start() + 1)}'.rjust(value_len, '.')[-value_len:]
                yield matches.group(1), new_id
    yield (ori_seq_ct, dig_seq_ct)

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f'author:\t{__author__}\nmail:\t{__mail__}\ndate:\t{__date__}\nversion:\t{__version__}'
    )
    parser.add_argument('-i', help='fa/fq文件或者数据库文件', dest='input', type=str, required=True)
    parser.add_argument('-e', help='酶切位点，可选 BcgI/CjePI/BsaXI/CspCI，默认为 BcgI', 
                      dest='enzyme', type=str, default='BcgI', choices=enzyme_pattern_dic.keys())  # 增加参数校验
    parser.add_argument('-o', help='输出文件前缀，建议使用样本名或者物种名', dest='output', type=str, required=False)
    parser.add_argument('-g', help='基因组文件，如果输入的是 CRAM 格式，该参数必须提供', dest='genome', type=str, required=False)
    parser.add_argument('-of', help='输出文件的压缩格式 marisa/gzip，默认为 marisa', 
                      dest='outfmt', type=str, choices=['marisa', 'gzip'], default='marisa')
    parser.add_argument('-l', help='marisa 值（id）长度，默认 50，建议不要更改', 
                      dest='value_len', type=int, default=50)
    parser.add_argument('--dump', help='将 marisa 转换成 fasta 格式（标准输出）', dest='dump', action='store_true')
    args = parser.parse_args()

    fmt = f'{args.value_len}c'
    if not args.output and not args.dump:
        report('ERROR', '缺少参数：给定 -o 将结果输出到文件，或者给定 --dump 将 marisa 转换成 fasta')

    if args.dump:
        trie = marisa_trie.RecordTrie(fmt).mmap(args.input)
        for tag in set(trie.keys()):
            for id_bytes in trie[tag]:
                # 修复bytes转字符串的方式
                id_str = ''.join([bytes([b]).decode('utf-8') for b in id_bytes]).lstrip('.')
                sys.stdout.write(f'>{id_str}\n{tag}\n')
    else:
        # 校验酶切位点是否存在（虽然已通过choices限制，但做双重保障）
        if args.enzyme not in enzyme_pattern_dic:
            report("ERROR", f"无效的酶切位点：{args.enzyme}，可选值：{list(enzyme_pattern_dic.keys())}")
        
        enzyme_pattern = re.compile(enzyme_pattern_dic[args.enzyme][1])
        
        # 收集酶切结果和统计信息
        results = extraction(args.input, args.value_len, enzyme_pattern, args.genome)
        marisa_key_lst, marisa_value_lst = [], []
        ori_seq_ct, dig_seq_ct = 0, 0
        
        for item in results:
            if isinstance(item, tuple) and len(item) == 2 and all(isinstance(x, int) for x in item):
                ori_seq_ct, dig_seq_ct = item
            else:
                seq, id = item
                if args.outfmt == 'marisa':
                    marisa_key_lst.append(seq)
                    marisa_value_lst.append(id.encode('utf-8'))
        
        if args.outfmt == 'marisa':
            trie = marisa_trie.BytesTrie(zip(marisa_key_lst, marisa_value_lst))
            trie.save(f'{args.output}.{args.enzyme}.fa.marisa')
            report("INFO", f"Marisa索引构建完成，共处理{ori_seq_ct}条序列，得到{dig_seq_ct}个酶切片段")
        else:
            with gzip.open(f'{args.output}.{args.enzyme}.fa.gz', 'wt') as OUT:
                results = extraction(args.input, args.value_len, enzyme_pattern, args.genome)
                for item in results:
                    if not isinstance(item, tuple) or len(item) != 2 or not all(isinstance(x, int) for x in item):
                        seq, id = item
                        OUT.write(f'>{id.lstrip(".")}\n{seq}\n')
            report("INFO", f"Gzip文件生成完成，共处理{ori_seq_ct}条序列，得到{dig_seq_ct}个酶切片段")
        
        with open(f'{args.output}.{args.enzyme}.dige.stat.xls', 'w') as STAT:
            STAT.write('sample\tenzyme\tinput_sequence_num\tenzyme_reads_num\tpercent\n')
            if ori_seq_ct == 0:
                percent = 0.0
            else:
                percent = round((dig_seq_ct / ori_seq_ct), 2)
            STAT.write(f'{args.output}\t{args.enzyme}\t{ori_seq_ct}\t{dig_seq_ct}\t{percent}\n')


if __name__ == "__main__":
    main()
