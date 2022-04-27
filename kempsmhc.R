library(ggplot2)
library(dada2) #make SURE you are running most updated version of R because the rcpp package required by dada2 doesn't work in previous version
library(reshape2)

path <- "~/Desktop/newkemps"

fnFs <- sort(list.files(path, pattern="_R1_SS.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_SS.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

plotQualityProfile(fnFs[1:2]) #look at quality profiles for first 2 samples for just R1reads, mean quality score is always above 30, higher past the first 25bp. indicates good quality trimming
plotQualityProfile(fnRs[1:2]) #R2 reads also look good, nothing below quality score of 30
#you can use dada2 to filter and trim, but I've already trimmed my reads with sickle/scythe, just change the path below to indicate that they are filtered

list.files(path) #make sure path lists all your sample files

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <-sample.names
names(filtRs) <-sample.names

out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = c(9,8), truncLen=c(0,0),
                   maxN=0, maxEE=c(0.1,0.1), truncQ=2, rm.phix=TRUE,
                   compress=TRUE, multithread=FALSE)
plotQualityProfile(filtFs[1:2]) #look at quality profiles for first 2 samples for just R1reads, mean quality score is always above 30, higher past the first 25bp. indicates good quality trimming
plotQualityProfile(filtRs[1:2]) #R2 reads also look good, nothing below quality score of 30

errF <- learnErrors(filtFs, multithread=TRUE) #32056257 total bases in 200352 reads from 35 samples will be used for learning the error rates
errR <- learnErrors(filtRs, multithread=TRUE) #32094297 total bases in 200352 reads from 35 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)#plot error rates for each possible nucleotide transition based on quality score. error rates should drop with increased quality, and black line should be good fit to observed points
plotErrors(errR, nominalQ = TRUE)

dadaFs<-dada(filtFs,err=errF, multithread = TRUE)
dadaRs<-dada(filtRs,err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <-mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose = TRUE)
head(mergers[[1]])

seqtab<-makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))

seqtab2<-seqtab[,nchar(colnames(seqtab)) %in% 190]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #leaves all 35 samples but only 310 variants remain, which is still a lot. if MOST were removed, might need to change upstream filters
sum(seqtab.nochim)/sum(seqtab2) 
write.table(seqtab.nochim, "newkemps_unfilt_variants.csv",sep=",",row.names = TRUE, col.names = TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

mod_seqtable<-seqtab.nochim
library("tidyverse")

mod_seqtable<-as.data.frame(t(as.matrix(seqtab.nochim))) #flip the dataframe so samples are cols and rows are variants, makes it easier to work with in dplyr
#as_tibble(mod_seqtable)
mod_seqtable$sums<- apply(mod_seqtable,1, function(i) sum(i>0)) #add column to dataframe showing how many samples have this variant
mod_seqtable<-dplyr::filter(mod_seqtable,sums>1) #keep only variants that are present in more than 1 sample, leaves 171 variants
mod_seqtable$total_calls<-apply(mod_seqtable,1, function(i) sum(i))
#mod_seqtable
#reads_by_samp<-as.data.frame(sum_reads_by_samp[c(1:35)],sample.names)

as.tibble(seqtab.nochim)
mod_seqtable<-as.data.frame(t(as.matrix(seqtab.nochim))) #flip the dataframe so samples are cols and rows are variants, makes it easier to work with in dplyr
as_tibble(mod_seqtable)

sum_reads_by_samp<-apply(mod_seqtable,2,function(i) sum(i)) #count total number of reads per sample
mod_seqtable$sums<- apply(mod_seqtable,1, function(i) sum(i>0)) #add column to dataframe showing how many samples have this variant

mod_seqtable$sums
mod_seqtable<-dplyr::filter(mod_seqtable,sums>1) #keep only variants that are present in more than 1 sample, leaves 171 variants
dim(mod_seqtable) #there actually weren't any variants in only one sample, so still have 181 variants
mod_seqtable$total_calls<-apply(mod_seqtable,1, function(i) sum(i))

reads_by_samp<-as.data.frame(sum_reads_by_samp[c(1:128)]) #put reads per sample in a dataframe with sample names

failed_samples<-ifelse (reads_by_samp<1000,"fail","pass") #consider failed samples as those than have fewer than 1k reads
failed_samples<-as.data.frame(failed_samples,stringsAsFactors = FALSE) #convert from character class to dataframe for filtering
failed_samples<-dplyr::filter(failed_samples, grepl("fail",`sum_reads_by_samp[c(1:128)]`)) #get list of failed samples

failed_samp_names<-row.names(failed_samples)  
#remove failed samples
mod_seqtable<-mod_seqtable[, !(names(mod_seqtable) %in% failed_samp_names)]
mod_seqtable<-mod_seqtable[, names(mod_seqtable) != "MHC_Sample72_lk_MA"]

mod_seqtable[]<- lapply(mod_seqtable,function(x) ifelse(x<100, 0, x)) #change calls to 0 if less than 100 supporting an allele
mod_seqtable$sums<- apply(mod_seqtable,1, function(i) sum(i>0)) #add column to dataframe showing how many samples have this variant
mod_seqtable<-dplyr::filter(mod_seqtable,sums>2) #same as above, redoin

seq_matrix<-as.matrix(mod_seqtable)#convert df to matrix to use prop table below
freq_matrix<-prop.table(seq_matrix,2) #for each variant in a sample, get the proportion of reads supporting that variant from total reads per individual
#replace all values with less than 0.03 proportion with 0, then recalculate sums for each variant, then filter variants
freq_matrix<-as.data.frame(freq_matrix)#convert back to dataframe
freq_matrix[]<-lapply(freq_matrix, function(x) ifelse(x<0.05, 0, x))
freq_matrix<-freq_matrix[c(1:121)]#remove last 3 columns that had sums/old frequency calcs
freq_matrix$sums<- apply(freq_matrix,1, function(i) sum(i>0))#add column with sums of number of samples per variant
freq_matrix<-dplyr::filter(freq_matrix,sums>1) #only keep variants in table if present in at least one individual
#this results in 60 variants 

final_var_list<-row.names(freq_matrix)
mod_seqtable<-mod_seqtable %>% rownames_to_column("row_names")
mod_seqtable2<-mod_seqtable[(mod_seqtable$row_names %in% final_var_list),]
row.names(mod_seqtable2)<-mod_seqtable2$row_names #convert row names back
mod_seqtable2<-mod_seqtable2[c(2:122)]#get rid of old sums and old row name column
mod_seqtable2$total_calls<-apply(mod_seqtable2,1, function(i) sum(i))
mod_seqtable2$sums<- apply(mod_seqtable2,1, function(i) sum(i>0))#add column with sums of number of samples per variant, REMEMBER that this column has one extra count in it because it is including the total calls column in count

total_reads<-sum(mod_seqtable$total_calls) #101231 total reads
mod_seqtable$variant_freq<-mod_seqtable$total_calls/total_reads*100

ggplot(mod_seqtable,aes(variant_freq))+
  geom_density()

no_alleles<-apply(mod_seqtable2,2,function(x) sum(x>0)) #get number of alleles per individual
no_alleles<-no_alleles[c(1:121)] #only keep allele numbers for first 34 columns which are the samples
no_alleles<-as.data.frame(no_alleles)
min(no_alleles$no_alleles) #1
max(no_alleles$no_alleles) #26
mean(no_alleles$no_alleles)#13

ggplot(no_alleles,aes(no_alleles))+
  geom_bar()+
  xlab('Number of sequence variants per individual')+
  ylab('Count Frequency')+
  theme_minimal()




