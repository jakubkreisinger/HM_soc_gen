library(dada2)
library(ggplot2)
library(ShortRead)

setwd("/media/kreising/DATA/data/RUN_16_3_2019_areny/02A.DEMULTI.16s/")

LIST<-list.files()
LIST
F_reads<-LIST[grep("pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("pair2.fastq.gz",LIST)]

system("zcat *pair1.fastq.gz > MERGED.fastq")
qp.f<-plotQualityProfile("MERGED.fastq")+ggtitle("Forward reads")
system("zcat *pair2.fastq.gz > MERGED.fastq")
qp.r<-plotQualityProfile("MERGED.fastq")+ggtitle("Rewerse reads")
system("rm MERGED.fastq")
qp.f
qp.r

ggsave(qp.f,filename = "/media/kreising/DATA/data/RUN_16_3_2019_areny/qp.f.pdf")
ggsave(qp.r,filename = "/media/kreising/DATA/data/RUN_16_3_2019_areny/qp.r.pdf")

sample.names<-gsub("-trimmed-pair1.fastq.gz","",F_reads)
sample.names<-gsub("assigned-","",sample.names)
tail(sample.names)

#POD TEMITO JMENY SE BUDOU UKLADAT FILTROVANE SOUBORY
filtFs <- paste0(sample.names, "_READ1_filt.fastq.gz")
filtRs <- paste0(sample.names, "_READ2_filt.fastq.gz")
tail(filtRs)
#FILTROVANI
#MAXIMALNI PREDPIKLADANY POCET CHYB - 1
for(x in 1:length(F_reads)) {
  print(sample.names[x])
  fastqPairedFilter(c(F_reads[x], R_reads[x]), c(filtFs[x], filtRs[x]),
                    maxN=0, maxEE=1, minQ=2,truncQ=2,
                    compress=TRUE, verbose=TRUE,
                    minLen = c(260,200),truncLen = c(260,200))
}


######COPY#####


FROM<-"/media/kreising/DATA/data/RUN_16_3_2019_areny/02A.DEMULTI.16s/"
TO<-"/media/kreising/DATA/data/RUN_16_3_2019_areny/FOR_DADA/"

SAMPLES<-list.files(FROM)
SAMPLES<-SAMPLES[grep("filt",SAMPLES)]
SAMPLES_new<-gsub("trus_READ","trus_Bendova02_READ",SAMPLES)

file.copy(from=paste0(FROM,SAMPLES),to=paste0(TO,SAMPLES_new))

#######################
#dada2 denoising#######
#######################

#These commands are executed on remote server (ca. 60G RAM, 400 CPU*hours)
library(dada2)

list.files()
fns <- list.files()
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("_READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("_READ2_filt.fastq.gz", fastqs)] 

sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#dereplication
derepFs <- derepFastq(fnFs,n = 1e+05, verbose=T)
derepRs <- derepFastq(fnRs,n = 1e+05, verbose=T)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#denoising
dadaFs <- dada(derepFs, multithread = 8,selfConsist = TRUE,MAX_CONSIST=25, pool = "pseudo")
dadaRs <- dada(derepRs, multithread = 8,selfConsist = TRUE,MAX_CONSIST=25, pool = "pseudo")

#merge denoised files
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10,maxMismatch=1,justConcatenate=F)

#abundance matrix
seqtab <- makeSequenceTable(mergers)

#save results
save(seqtab,file="/storage/brno2/home/kreising/otutab_16s_pseudo.R")
write.table(seqtab,"/storage/brno2/home/kreising/otutab_pseudo.txt",sep="\t")








