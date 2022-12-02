# MAP2B

----------------------------
**M**et**A**genomic **P**rofiler based on type **IIB** restriction site

### Why do I need to run MAP2B?
Accurate species identification and abundance estimation are critical for the interpretation of whole metagenome sequencing (WMS) data. Numerous computational methods, broadly referred to as metagenomic profilers, have been developed to identify species in microbiome samples by classification of sequencing reads and quantification of their relative abundances. Yet, existing metagenomic profilers typically suffer from false-positive identifications and consequently biased relative abundance estimation. Indeed, false positives can be accounted for more than 90% of total identified species. Here, we present a new metagenomic profiler MAP2B to resolve those issues. 

For more details, please see [Eliminate false positives in metagenomic profiling based on type IIB restriction sites](https://www.biorxiv.org/content/10.1101/2022.10.24.513546v1).

## The workflow of MAP2B
Instead of directly estimating the relative abundances of the species through aligning reads against the whole microbial genome or marker genes as existing metagenomic profilers do, we use the following two-round reads alignment strategy:

(**A**) For any input WMS data, 2b tags can be extracted by 2B in silico digestion.   
(**B**) WMS-originated 2b tags will be mapped against a preconstructed unique 2b tag database, which contains ~50,000 identifiable species.  
(**C**) In the 1st round of reads alignment, genome coverage, taxonomic count, sequence count and G score are calculated for each species.  
(**D**) The four features above will be passed into a preconstructed false positive recognition model.  
(**E**) A high-precision species identification result will be generated.  
(**F**) A sample-dependent unique 2b tag database will be constructed based on the species identification result.  
(**G**) In the second round of reads alignment, we estimate [taxonomic abundance](https://www.nature.com/articles/s41592-021-01141-3) for each species.  

 ![workflow](MAP2B_workflow.png)
 
## Installation
 
 ### System requirements
 
 #### Dependencies
All scripts in MAP2B are programmed by Perl and Python, and execution of MAP2B is recommended in a conda environment. This program could work properly in the Unix systems, or Mac OSX, as all required packages can be appropreiately download and installed.  
 #### Memory usage
**> 120G RAM** is required to run this pipeline.  
 ### Download the pipeline
 * Clone the latest version from GitHub (recommended):  
 
   `git clone https://github.com/sunzhengCDNM/MAP2B/`  
   `cd MAP2B`
   
    This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added.
 * Alternatively, directly download the whole GitHub repo without installing GitHub:
 
   `wget https://github.com/sunzhengCDNM/MAP2B/archive/master.zip`  
   `unzip master.zip`  
   `cd MAP2B-master`
   
 ### Install MAP2B in a conda environment 
 * Conda installation  
   [Miniconda](https://docs.conda.io/en/latest/miniconda.html) provides the conda environment and package manager, and is the recommended way to install MAP2B. 
 * Create a conda environment for MAP2B pipeline:  
   After installing Miniconda and opening a new terminal, make sure you’re running the latest version of conda:
   
   `conda update conda`
   
   Once you have conda installed, create a conda environment with the yml file `tools/MAP2B-20221128-conda.yml`.
   
   `conda env create -n MAP2B-2022 --file tools/MAP2B-20221128-conda.yml`
   
 * Activate the MAP2B conda environment by running the following command:
 
   `source activate MAP2B-2022` or `conda activate MAP2B-2022`
   
   Make sure the conda environment of MAP2B has been activated by running the above command before you run MAP2B everytime.  

 * Run the following command to download the preconstructed unique 2b tag database:
 
   `perl scripts/Download_MAP2B_DB_GTDB.pl database`
   
   Resuming from a breakpoint is supported during database download.  
   Now, everything is ready for MAP2B :), Let's get started.
 
## Using MAP2B
 
### Quick start
MAP2B is a highly automatic pipeline, and only a few parameters are required for the pipeline.
* We prepared a real pair-end sequencing data of a MOCK community:  
 
   `cd example`  
   `mkdir -p data/`  
   `wget -t 0 -O data/shotgun_MSA-1002_1.fq.gz https://figshare.com/ndownloader/files/38346149/shotgun_MSA-1002_1.fq.gz`  
   `wget -t 0 -O data/shotgun_MSA-1002_2.fq.gz https://figshare.com/ndownloader/files/38346155/shotgun_MSA-1002_2.fq.gz`  
 
* After downloading the sequencing data, we can finally run MAP2B:  
 
   `python3 ../bin/MAP2B.py -i data.list`

    In `data.list` you can learn how to prepare your input data, both single-end and paired-end data can be used as input.  
    
```
sample1<tab>shotgun1_left.fastq(.gz)<tab>shotgun1_right.fastq(.gz)
sample2<tab>shotgun2_left.fastq(.gz)<tab>shotgun2_right.fastq(.gz)
sample3 ...

```

### Parameters
The main program is `bin/MAP2B.py` in this repo. You can check out the usage by printing the help information via `python3 bin/MAP2B.py -h`.

```
usage: MAP2B.py [-h] -i INPUT [-e ENZYME] [-d DATABASE] [-p PROCESSES]
                [-o OUTPUT]

optional arguments:
  -h, --help    show this help message and exit
  -i INPUT      The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. like 
                  sample <tab> shotgun.1.fq(.gz) (<tab> shotgun.2.fq.gz)
  -e ENZYME     enzyme, default 13 for CjePI, choose from
                  [1]CspCI  [5]BcgI  [9]BplI     [13]CjePI  [17]AllEnzyme
                  [2]AloI   [6]CjeI  [10]FalI    [14]Hin4I
                  [3]BsaXI  [7]PpiI  [11]Bsp24I  [15]AlfI
                  [4]BaeI   [8]PsrI  [12]HaeIV   [16]BslFI
  -d DATABASE   Database path for MAP2B pipeline, default MAP2B/database
  -p PROCESSES  Number of processes, note that more threads may require more memory, default 1
  -o OUTPUT     Output directory, ./MAP2B_result

author: Liu Jiang, Zheng Sun
mail: jiang.liu@oebiotech.com, spzsu@channing.harvard.edu
Last update: 2022/11/22 11:21:47
version:  1.0.0
```

## Reference
 * Sun, Z., Liu, J., Zhang, M. et al. Eliminate false positives in metagenomic profiling based on type IIB restriction sites. bioRxiv 2022.10.24.513546. https://doi.org/10.1101/2022.10.24.513546  
 * Sun, Z., Huang, S., Zhu, P. et al. Species-resolved sequencing of low-biomass or degraded microbiomes using 2bRAD-M. Genome Biol 23, 36 (2022). https://doi.org/10.1186/s13059-021-02576-9  
 * Sun, Z., Huang, S., Zhang, M. et al. Challenges in benchmarking metagenomic profilers. Nat Methods 18, 618–626 (2021). https://doi.org/10.1038/s41592-021-01141-3  
 * Wang, S., Liu, P., Lv, J. et al. Serial sequencing of isolength RAD tags for cost-efficient genome-wide profiling of genetic and epigenetic variations. Nat Protoc 11, 2189–2200 (2016). https://doi.org/10.1038/nprot.2016.133  
 * Wang, S., Meyer, E., McKay, J. et al. 2b-RAD: a simple and flexible method for genome-wide genotyping. Nat Methods 9, 808–810 (2012). https://doi.org/10.1038/nmeth.2023  

## Acknowledgement
This work was supported by the National Institutes of Health grant number R01AI141529, R01HD093761, RF1AG067744, UH3OD023268, U19AI095219, U01HL089856, and the Charles A. King Trust Postdoctoral Fellowship. 
