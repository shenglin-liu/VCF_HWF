#!/usr/bin/env Rscript

# Hardy-Weinberg filtering of a VCF file.

# Usage (bash): [path_to]/VCF_HWF.r $vcf_file $popmap_file $min_pop_HWE $out_vcf

# Requirements:
#	genotype fields start with GT tag.
#	biallelic, 0 representing REF and 1 for ALT.

input<-commandArgs(trailingOnly=TRUE)
f.vcf<-input[1]
f.popmap<-input[2]
minPopHWE<-input[3]
f.out<-input[4]

## Loading the genotyping data in substrings
vcf<-read.table(f.vcf,sep="\t",stringsAsFactors=FALSE)
vcf.gt<-as.matrix(vcf[,-c(1:9)])
nrow.vcf<-nrow(vcf.gt)
ncol.vcf<-ncol(vcf.gt)
chrom1<-matrix(as.integer(substring(vcf.gt,1,1)),nrow.vcf,ncol.vcf) #Extracts a substring (0/0 - it takes the first 0)
chrom2<-matrix(as.integer(substring(vcf.gt,3,3)),nrow.vcf,ncol.vcf) #Extracts a substring (0/0 - it takes the second 0)
chrom<-chrom1+chrom2 #Creates a matrix containing only information about Homoz (0 or 2) or Hetz(1)

## Loading popmap
popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=FALSE)[,2]
pops<-unique(popmap)
n<-length(pops)

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
FissqrtN<-Fis*sqrtN
Boolean<-(FissqrtN>1.96)|(FissqrtN<(-1.96)) ##Makes a true/false matrix
Blacklist<-which(rowSums(Boolean,na.rm=TRUE)>(n-minPopHWE)) # The SNPs where more than (n-minPopHWE) populations are not in HW equilibrium
write.table(vcf[-Blacklist,],f.out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE) #now we have a HW-filtered vcf file
