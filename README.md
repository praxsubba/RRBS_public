RRBS-Workflow

Mapping RRBS Reads to the Reference Genome

bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam -D C-CGG -D T-CGA -w 100 -v 0.08 -r 0 -p 4 -n 0 -s 12 -S 0 -f 5 -q 0 -u -V 2.

#Reference genome: GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff

Minimap2

#This is used to map cDNA sequences from Dong et al., 2009 to the zebrafinch transcriptome GCF_008822105.2_bTaeGut2.pat.W.v2_rna.fna

minimap2 -I 13G /scratch1/psubba/GCF_008822105.2_bTaeGut2.pat.W.v2_rna.fna.gz /scratch1/psubba/sb_array_s eq.FASTA > approx-mapping.txt

How to get filter apprx-mapping.txt with significant transcript_IDs from fdr clone of pairwise comparisons

setwd("~/Box/Prakrit Subba rotation George Lab/microarray_experiment")

#Make sure that if there are two words, you delete the space to make it one word. Example, Histone H3 should be Histone.H3.

a <- read.table(file = "AASA_gene_list.txt")

b <- read.table(file = "approx-mapping_final.txt")

install.packages("dplyr")

library("dplyr")

colnames(a) <- "transcript.id"

colnames(b) <- c("transcript.id","tavish")

c <- dplyr::inner_join(a,b,by="transcript.id")

write.table(c,file="AASA_sig_transcript_ID.txt", sep = "\t",row.names = FALSE)

Replace AASA with whatever the comparison might be

Manually remove the first column line before transferring to cluster

#Remove double quotes

sed 's/"//g' AASA_sig_transcript_ID.txt > AASA_sig_transcript_ID_2.txt

If you used the R script with dplyr, then skip to the lines after "This is the second file that will be merged"

Fixing Gene IDs / methratio Python Error

cut -f 1,6 AASS_sig_transcript_ID.txt >AASS_sig_merge_1.txt

#AASS_sig_transcript has columns 1 and 6 from the output of minimap2

#AASS_sig_merge_1.txt is the first file to be used for the merge step

#This is the second file that will be merged

If you used R/dplyr, start below this. Also note that if you've obtained "gene_tx_id_merge_2_final.txt" once, you can resuse it

grep ">" /scratch1/psubba/GCF_008822105.2_bTaeGut2.pat.W.v2_rna.fna > headers.txt #extract fasta headers

cat headers.txt | cut -d " " -f1 |sed 's/^>//' > tx_id.txt

cat headers.txt | rev | cut -d ")" -f2 | cut -d "(" -f1 | rev > gene_id.txt

#created the second file to be merged with two columns (genes and transcript IDs)

paste gene_id.txt tx_id.txt | column -s $'\t' -t > gene_tx_id_merge_2_final.txt

#In the directory /home/psubba/fixing_gene_id, files were merged to have: SB Clone ID, Gene name, Transcript ID.

Merging AASS_sig_merge_1.txt with Final_merged.txt

##Script with significant genes: AASS_sig_merge_1.txt

##Script with gene IDs and transcript IDs: gene_tx_id_merge_2_final.txt

awk '{print $2}' AASS_sig_merge_1.txt | while read line: do grep $line gene_tx_id_merge_2_final.txt; done > Final_merged.txt

Alternatively use the following R script for merging

setwd("~/Box/Prakrit Subba rotation George Lab/microarray_experiment")

a <- read.table(file = "AASS_sig_transcript_ID_2.txt")

b <- read.table(file = "gene_tx_id_merge_2_final.txt")

install.packages("dplyr")

library("dplyr")

colnames(a) <- c("1","2")

colnames(b) <- c("1","2")

c <- dplyr::left_join(a,b,by="2")

colnames(c) <- c("clone_ID","transcript_ID","gene")

write.table(c,file="Final_merged_AASS.txt",sep = "\t",row.names = FALSE)

If you used the R script, then, skip to finding unique genes

awk '{print $1}' AASS_sig_merge_1.txt > first.txt

pr -tmJ Final_merged.txt first.txt > output.txt

#Output.txt contains: SB Clone ID, Gene name, Transcript ID

#Find unique (i.e., non-repeating/deduplicated) genes in output.txt

If you used the R script, pick up below this.

sed 's/"//g' Final_merged_AASA.txt > Final_merged_AASA_2.txt

awk '{print $3}' Final_merged_AASA_2.txt | sort | uniq -u > unique_genes_AASA.txt #Use this if you USED the R script

awk '{print $1}' output.txt | sort | uniq -u > unique_genes.txt #Use this if you did not use the R script

#unique_genes.txt has genes that are significant and this will be used to filter the RDS file from Rtracklayer. Why? Because R track layer has all the genes and the information from that. We are just interested in these genes in unique_genes.txt

Rtracklayer Script for getting gene intervals (as well as 10,000 bp before start site or upstream of every gene) title: "extracting gene intervals from GFF file, zebra finch bTaeGut1_v1_p" output: html_notebook

setwd("~/Box/Prakrit Subba rotation George Lab/rtracklayer_corrected")
library("rtracklayer")
library("tidyverse")
readGFF("GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff")%>%head()
I changed this code to select "gene" instead of "mRNA", "lncRNA", "transcript" Now there should be only one line for each gene (instead of the multiple transcripts)

my_tags <- c("Name", "Dbxref","gene")
my_columns <- c("seqid", "start", "end", "strand", "type")
my_filter<-list(type="gene")
dat<-readGFF("GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff",tags=my_tags,columns=my_columns,filter=my_filter)
head(dat)
I made a new column for interval_start and interval_stop, so I could increase the upstream interval to capture potential regulatory sites

dat$interval_start<-dat$start
dat$interval_stop<-dat$end
dat
Add 10kb upstream of the start site gene intervals always have start and stop listed smallest to largest (like in a bed file) for genes on the positive strand, start = start-10000 for genes on the negative strand, stop = stop + 10000

as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-10000,start))%>%  #if positive strand, subtract 10000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 10kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+10000,end))%>%   #if negative strand, add 10000 to the end
  select(seqid,interval_start,interval_stop)
save the table as a txt file

as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-10000,start))%>%  #if positive strand, subtract 10000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 10kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+10000,end))%>%   #if negative strand, add 10000 to the end
  select(seqid,interval_start,interval_stop)%>%write.table("genewise_bed_bTaeGut_v1_p_v2.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
export a version for mapping methylation sites back to genes (all columns)

#print to screen
as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-10000,start))%>%  #if positive strand, subtract 10000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 10kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+10000,end))

#write to file
as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-10000,start))%>%  #if positive strand, subtract 10000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 10kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+10000,end))%>%
  saveRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS")
test <- readRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS") #This reads the RDS file into a data frame test
unique_genes_AASA <- read.table(file = "unique_genes_AASA.txt") #This reads the txt file in a data frame unique_genes
colnames(unique_genes_AASA) <- "gene" #This step changes the column name/header to "gene" for the unique_genes data frame
TEST_2 <- dplyr::semi_join(test, unique_genes_AASA, by="gene") #This step matches filters the RDS for only the unique or significant genes in unique_genes
saveRDS(TEST_2, file="AASA_bTaeGut_v2p_gene_intervals_plus_10K_final.RDS") #saves the RDS file
write.table(test,file="AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2.txt",sep = "\t",
            row.names = FALSE) #this txt file can be used as a BED file in the following steps
sessionInfo()
Preparing BED files to extract alignments from BSMAP output files

cut -f 1,9,10 AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2.txt > AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut.txt #extracted the newid, gene interval start and stop sites

tail -n +2 AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut.txt > AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_2.txt #removed the first row from the file

cat AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_2.txt | sed 's/"//g' > AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_3.bed #removed double quotes from the file

OR

cat AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_2.txt | sed 's/"//g' > AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_3.txt #removed double quotes from the file

Extracting alignments from BSMAP output files using samtools view

samtools view -b -L /home/psubba/fixing_gene_id/AASA_bTaeGut_v2p_gene_intervals_plus_10K_v2_cut_3.bed /home/psubba/BSMAP_output_BAM/$_output.bam > $_interval.bam

Running methratio.py from BSMAPz

#PBS -N Methratio_G1FA211 #PBS -l select=1:ncpus=12:mem=62gb:interconnect=fdr,walltime=24:00:00 #PBS -m abe #PBS -j oe

cd $PBS_O_WORKDIR python ~/miniconda2/envs/methratio/bin/methratio.py -o $_methratio.txt -d /scratch1/psubba/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna -z -x CG /home/psubba/BSMAP_output_BAM/$_interval.bam

MethylKit

title: "RRBS targeted methyl kit Final by Prakrit" output: html_notebook

#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")
library("methylKit")
setwd("~/Box/Prakrit Subba rotation George Lab/methratio_files_final/")
file.list = list("G1FA160_methratio.txt", "G1FA211_methratio.txt", "G1FA237_methratio.txt", "G1SI152_methratio.txt", "G1SI199_methratio.txt", "G1SI246_methratio.txt", "G2FA149_methratio.txt", "G2FA209_methratio.txt", "G2FA222_methratio.txt", "G2SI146_methratio.txt", "G2SI188_methratio.txt", "G2SI218_methratio.txt")
myobj=methRead( file.list,pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2, coverage.col=6,strand.col=3,freqC.col=5 ),
                sample.id=list("G1FA160", "G1FA211", "G1FA237", "G1SI152", "G1SI199", "G1SI246", "G2FA149", "G2FA209", "G2FA222", "G2SI146", "G2SI188", "G2SI218"),assembly="taeGut1",treatment=c(0,0,0,1,1,1,0,0,0,1,1,1))
#get methylation statistics

for (i in 1:12) {
  getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
}
#filter by coverage and get coverage stats

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=NULL)

for (i in 1:12){
  getCoverageStats(filtered.myobj[[i]],plot=TRUE,both.strands=FALSE)
}
meth=unite(filtered.myobj, destrand=FALSE)
meth
meth_destrand=unite(filtered.myobj, destrand=TRUE) 
##see that destrand=TRUE increases coverage
meth_destrand
#see how samples cluster

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE,sd.threshold = .90)
clusterSamples(meth_destrand, dist="correlation", method="ward", plot=TRUE,sd.threshold = .90)
#plot Principal Components

PCASamples(meth,adj.lim = c(.4, 1), sd.threshold = .90) #not sure about how to set this value?
PCASamples(meth_destrand,adj.lim = c(.5, 1), sd.threshold = .90) #not sure about how to set this value?

#get a methylDiff object containing the differential methylation statistics and locations for regions or bases
myDiff=calculateDiffMeth(meth)
#write.csv(getData(myDiff),"summary.csv")
getData(myDiff)
myDiff_d=calculateDiffMeth(meth_destrand)
getData(myDiff_d)
# get all differentially methylated bases, difference greater than 25
myDiff10p=getMethylDiff(myDiff,qvalue=.01)
getData(myDiff10p)


myDiff10p_destranded=getMethylDiff(myDiff_d,qvalue=.01)
getData(myDiff10p_destranded)







sessionInfo()
