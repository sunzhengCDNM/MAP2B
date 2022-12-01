#!/usr/bin/perl
#Author:zhangrongchao, zhangrongchaoxx@163.com
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

$ARGV[0] ||="MAP2B_DB_GTDB";

my $outdir=$ARGV[0];#下载目录

$outdir=abs_path($outdir);
&CheckDir("$outdir");


my @a=('abfh_classify', 'CjePI.species.stat', 'CjePI.species.uniq.stat');#分类表及统计文件
my @b=('CjePI.species', 'CjePI.species.uniq');#需要下载的库文件

my %hash_path=(
	'abfh_classify'          =>['https://figshare.com/ndownloader/files/37249060/abfh_classify_with_speciename.txt.gz',],
	'CjePI.species'          =>['https://figshare.com/ndownloader/files/37252444/CjePI.species.fa.gz0',
	                            'https://figshare.com/ndownloader/files/37248751/CjePI.species.fa.gz1',
	                            'https://figshare.com/ndownloader/files/37248892/CjePI.species.fa.gz2',
	                            'https://figshare.com/ndownloader/files/37249057/CjePI.species.fa.gz3',
	                            'https://figshare.com/ndownloader/files/37259224/CjePI.species.fa.gz4',
	                            'https://figshare.com/ndownloader/files/37252729/CjePI.species.fa.gz5',
	                            'https://figshare.com/ndownloader/files/37252894/CjePI.species.fa.gz6',],
	'CjePI.species.stat'     =>['https://figshare.com/ndownloader/files/38340911/CjePI.species.stat.xls'],
	'CjePI.species.uniq'     =>['https://figshare.com/ndownloader/files/38334476/CjePI.species.uniq.fa.gz1',
	                            'https://figshare.com/ndownloader/files/38334494/CjePI.species.uniq.fa.gz2',],
	'CjePI.species.uniq.stat'=>['https://figshare.com/ndownloader/files/38334467/CjePI.species.uniq.stat.xls',],
	);

my %hash_md5=(
	'abfh_classify'          =>['6f501521429b7f0fd6ec3a33ffb5b7c8',],
	'CjePI.species'			 =>['a847bce6f10631382aca8899dfc0bcf4',
	                            'e074cc01c0bcc6a7e41cabdcc29dd388',
	                            '837ad8f24281a853cd2f45eff0d6c3e0',
	                            'a5d6a1f9b6d27c2bc476e5c7ff73aa6a',
	                            'ce7b6214b8344dc33f53d8b735f9a02b',
	                            '0ba9bc286fddec34a3e5248fb9b8de10',
	                            '99431d50b2aab1eefe2ebcca873e1168',],
	'CjePI.species.stat'     =>['edfbdb0a61b77ee274518da263ee0a8f',],
	'CjePI.species.uniq'     =>['36a6c842c286de5d76933ba72854acb8',
	                            '6e9c44aa5fb6a0ec915e8c8fdd423b9d',],
	'CjePI.species.uniq.stat'=>['037324b8abab2dcafb36cdd48e24c1eb',],
	);

#合并后文件md5
my %complete_md5=(
	'CjePI.species'           =>'83065fce3b160a2ad73aadced0b3582c',
	'CjePI.species.uniq'      =>'15d789c6d0db56e4e30f820399b36163',
	);

#download abfh_classify
for my $i(@a){
	my @tmp=split /\//,$hash_path{$i}[0];
	my $url=join("/",@tmp[0..$#tmp-1]);
	my $name=$tmp[-1];
	my $file_md5;#下载的文件的MD5值
	while(1){
		if(-e "$outdir/$name"){
			chomp($file_md5=`md5sum $outdir/$name`);
			$file_md5=(split /\s+/,$file_md5)[0];
		}
		if(-e "$outdir/$name" && $file_md5 eq $hash_md5{$i}[0]){
			print STDOUT "File $name has been downloaded.\n";
			last;
		}else{
			`wget -t 0 -O $outdir/$name $url`;
		}
	}
}

#下载数据库文件
for my $i(@b){
	my $cat="";
	while(1){
		my $md5;
		if(-e "$outdir/$i.fa.gz"){#存在完成文件
			chomp($md5=`md5sum $outdir/$i.fa.gz`);
			$md5=(split /\s+/,$md5)[0];
		}
		if(-e "$outdir/$i.fa.gz" && $md5 eq $complete_md5{$i}){
			print STDOUT "File $i.fa.gz hash been downloaded.\n";
			`rm -rf $cat`;
			last;
		}else{
			for my $j(0..$#{$hash_path{$i}}){#循环每个文件
				my @tmp=split /\//,$hash_path{$i}[$j];
				my $url=join("/",@tmp[0..$#tmp-1]);
				my $name=$tmp[-1];
				my $file_md5;#下载的文件的MD5值
				while(1){
					if(-e "$outdir/$name"){
						chomp($file_md5=`md5sum $outdir/$name`);
						$file_md5=(split /\s+/,$file_md5)[0];
					}
					if(-e "$outdir/$name" && $file_md5 eq $hash_md5{$i}[$j]){
						print STDOUT "File $name has been downloaded.\n";
						$cat .=" $outdir/$name";
						last;
					}else{
						`wget -t 0 -O $outdir/$name $url`;
					}
				}
			}
			`cat $cat > $outdir/$i.fa.gz`;
		}
	}
}

print STDOUT "Congratulations! All databases have been downloaded.\n";

sub CheckDir{
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}
