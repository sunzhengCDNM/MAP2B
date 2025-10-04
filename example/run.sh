# conda activate MAP2B.1.5
# pre data
mkdir -p data/ && \
wget --user-agent="Mozilla/5.0" -t 3 -O data/shotgun_MSA-1002_1.fq.gz https://figshare.com/ndownloader/files/38346149/shotgun_MSA-1002_1.fq.gz && \
wget --user-agent="Mozilla/5.0" -t 3 -O data/shotgun_MSA-1002_2.fq.gz https://figshare.com/ndownloader/files/38346155/shotgun_MSA-1002_2.fq.gz && \
# run demo
time python3 ../bin/MAP2B.py -i data.list
