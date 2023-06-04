#!/usr/bin/env Rscript

# Hardy-Weinberg filtering of a VCF file.

# Usage (bash): [path_to]/VCF_HWF.r $vcf_file $popmap_file $min_pop_HWE $out_vcf
#	vcf_file: the input file in VCF format.
#	popmap_file: a population map file in the same format as in Stacks (https://catchenlab.life.illinois.edu/stacks/manual/).
#	min_pop_HWE: minimum number of populations in HWE for a SNP to be kept.
#	out_vcf: name of the new VCF file to output.

# Requirements:
#	genotype fields start with GT tag.
#	biallelic, 0 representing REF and 1 for ALT.

input<-commandArgs(trailingOnly=TRUE)
f.vcf<-input[1]
f.popmap<-input[2]
minPopHWE<-as.integer(input[3])
f.out<-input[4]

vcf_hwf<-function(vcf,popmap,minPopHWE)
{
	pops<-unique(popmap)
	n<-length(pops)

	## Parsing the genotyping data in substrings
	vcf.gt<-as.matrix(vcf[,-c(1:9),drop=FALSE])
	nrow.vcf<-nrow(vcf.gt)
	ncol.vcf<-ncol(vcf.gt)
	defaultW<-getOption("warn")
	options(warn=-1)
	chrom1<-matrix(as.integer(substring(vcf.gt,1,1)),nrow.vcf,ncol.vcf) #Extracts a substring (0/0 - it takes the first 0)
	chrom2<-matrix(as.integer(substring(vcf.gt,3,3)),nrow.vcf,ncol.vcf) #Extracts a substring (0/0 - it takes the second 0)
	options(warn=defaultW)
	chrom<-chrom1+chrom2 #Creates a matrix containing only information about Homoz (0 or 2) or Hetz(1)

	## Loop for observed and expected population-based heterozygosity
	Hobs<-Hexp<-Nsam<-matrix(0,nrow.vcf,n)
	colnames(Hobs)<-colnames(Hexp)<-colnames(Nsam)<-pops
	for(pop in pops)
	{
		index<-popmap==pop
		chrom.pop<-chrom[,index,drop=FALSE]
		Hobs[,pop]<-rowMeans(chrom.pop==1,na.rm=TRUE)
		af1<-rowMeans(chrom.pop/2,na.rm=TRUE)
		af0<-1-af1
		Hexp[,pop]<-2*af1*af0
		Nsam[,pop]<-rowSums(!is.na(chrom.pop))
	}

	## Hardy-Weinberg filtering
	sqrtN<-Nsam^0.5
	Fis<-(Hexp-Hobs)/Hexp
	Fis[Hexp==0]<-0 #monomorphic cases are set to 0 (otherwise, they will be NA)
	FissqrtN<-Fis*sqrtN
	Boolean<-(FissqrtN<=1.96)&(FissqrtN>=(-1.96)) ##Makes a true/false matrix
	Whitelist<-which(rowSums(Boolean,na.rm=TRUE)>=minPopHWE) #The SNPs where at least minPopHWE populations are in HWE
	vcf[Whitelist,]
}

## Loading popmap
popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=FALSE)[,2]

## Loading and processing the VCF chunk by chunk through connection (it can handle large VCF files)
c.vcf<-file(f.vcf,"r")
c.out<-file(f.out,"w")
chunk<-2e5 #Chunk size, i.e., number of lines to read in each cycle
repeat
{
	vcf<-try(read.table(c.vcf,sep="\t",stringsAsFactors=FALSE,colClasses="character",nrows=chunk),silent=TRUE)
	if(class(vcf)!="data.frame")break
	out<-vcf_hwf(vcf,popmap,minPopHWE)
	write.table(out,c.out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE) #now we have a HW-filtered vcf file
}
close(c.vcf)
close(c.out)

# Be careful:
#	one pop.
#	one line or one column cases.
#	difference between a blacklist of ">=(n-minPopHWE)" and a whitelist of ">=minPopHWE" due to NA.
#	sample size = 0 due to missing.
#	Fis = (0-0)/0 due to monomorphism.
# All checked; the script is robust.
