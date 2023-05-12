###########################################################################
# Function: Abundance table simulation and sequencing reads simulation
# Call: Rscript simulation_abd.R type[gut/oral/skin/buliding/ocean/vaginal] [sample#] [min_species#] [max_species#] [distribution_type] [enable_sequencing_reads_generation][T/F] [reads#]
# Example: Rscript simulation_abd.R gut 100 100 200 L T 100000
# Explanation: Generating an abudance table and sequencing data from zero-inflated log-normal distribution with 100 samples, species number is randomly selected range from 100 to 200 in each sample.
# Genome_ref and Wgsim are necessary for the script
# Species number included in database:
# ----------------------------------------------------------------------------------
#	Gut	Oral	Skin	Vaginal	Ocean	human	building
# ----------------------------------------------------------------------------------
# 	409	445	272	158	109	1287	349
# ----------------------------------------------------------------------------------
# For more information about GCF accesion number, please refers to GCF_annotation.
# L: lognormal distribution, R: random distribution, N: one-side normal distribution, E:even distribution.
# Last update: 2020-04-01, Zheng Sun 
###########################################################################

setwd('./')

#Args input
args=commandArgs(T)
if(is.na(args[1])) stop("Please input the parameters e.g. Rscript simulation_abd.R gut 100 100 200 L F")
if(is.na(args[5])) stop("Please input the parameters e.g. Rscript simulation_abd.R gut 100 100 200 L F")
if(is.na(args[6])) stop("Please input the parameters e.g. Rscript simulation_abd.R gut 100 100 200 L F")

args[1] #type
args[2] #sample#
args[3] #min species#
args[4] #max sample#
args[5] #distribution type
args[6] #if enable sequencing data generating
args[7] #reads#
outputname <- paste(args[2],"_",args[1],"_abd.txt",sep="") #prefix assignment
print(paste("Generating an abundance table named",outputname ,"from a specifed distribution with", args[2], args[1],"samples, species number is randomly selected range from",args[3] ,"to", args[4],"in each sample."))
outputfolder <- paste(args[2],"_",args[1],"_samples_simulation",sep="")  #output folder name
system(paste("mkdir ",outputfolder,sep="")) #output folder creation

#Args assignment
genome_length <- read.table(paste("Genome_ref/",args[1],sep=""),header = T,sep="\t")
sample_num <- args[2]
min_spp <- args[3] #assigning the minimum nubmer of spp in simulation data
max_spp <- args[4] #assigning the maximum number of spp in simulation data
max_species_db <- dim(genome_length)[1] #calculate the maximum number of spp in database
smo_taxonabd <- cbind() # original simulated matrix (taxonomic abd)
smo_seqabd <- cbind() # original simulated matrix (sequence abd)

if (args[5]=="L") { #log-normal distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(rlnorm(species_num),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args[5]=="N") { #one-side normal distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(abs(rnorm(species_num)),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args[5]=="R") { #random distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(runif(species_num),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args[5]=="E") { #even distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(rep(1, times=(species_num)),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else {
stop("Please assign a type of distribution for abundance simulation.")
}

smo_seqabd_origin <- smo_seqabd*(as.numeric(args[7])) #for generating sequencing reads

#generating tables
species_ID <- c("SpeciesID",as.character(genome_length[,3])) #factor to character
smo_taxonabd <- rbind(paste("sample",seq(from=1,to=args[2]),sep="_"),smo_taxonabd) #header added
smo_taxonabd <- cbind(species_ID,smo_taxonabd) #rownames added
smo_seqabd <- rbind(paste("sample",seq(from=1,to=args[2]),sep="_"),smo_seqabd ) #header added
smo_seqabd <- cbind(species_ID,smo_seqabd ) #rownames added

write.table(smo_taxonabd,file=paste(outputfolder,"/",args[2],"_",args[1],"_taxonomic_abd.txt",sep=""),quote=F,,row.names=F,col.names = F,sep = "\t")
write.table(smo_seqabd,file=paste(outputfolder,"/",args[2],"_",args[1],"_sequence_abd.txt",sep=""),quote=F,,row.names=F,col.names = F,sep = "\t")

if (args[6]=="T") { #if enable generating sequencing reads
seq_abd_for_wgsim <- as.data.frame(cbind(as.matrix(smo_seqabd_origin), as.character(genome_length[,4])))  #merge seq abd and its genome name
print(seq_abd_for_wgsim)
for (lines in 1:(dim(seq_abd_for_wgsim)[2]-1)) { # read seq abd table and split the table into samples line by line
	temp <- cbind(seq_abd_for_wgsim[lines],seq_abd_for_wgsim[dim(seq_abd_for_wgsim)[2]]) # creat an temp table to save the current line splited by seq table
	temp <- temp[which(rowSums(temp==0)==0),] # discard the lines with zero
	system(paste("mkdir ",outputfolder,"/sample",lines,sep="")) # creat an folder to save the simulation data
for (i in 1:dim(temp)[1]) { # in the current sample, [i,1] is the number of the reads, [i,2] is the name of genome
	system(paste("Wgsim/wgsim -N ",ceiling(as.numeric(as.character((temp[i,1]))))," -1 150 -2 150 Genome_ref/",temp[i,2]," ",outputfolder,"/sample",lines,"/",temp[i,2],".Left.reads "," ",outputfolder,"/sample",lines,"/",temp[i,2],".Right.reads",sep=""))
	# Wgsim is used to generate the simulation data
}
	system(paste("cat ",outputfolder,"/sample",lines,"/*.Left.reads > ",outputfolder,"/sample",lines,".left.fq",sep=""))
	system(paste("cat ",outputfolder,"/sample",lines,"/*.Right.reads > ",outputfolder,"/sample",lines,".right.fq",sep=""))
	system(paste("rm -rf ",outputfolder,"/sample",lines,sep=""))
} 
} else {
print("No sequencing reads were simulated.")
}


