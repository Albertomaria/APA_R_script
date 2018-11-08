library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(tidyverse)
library(TransView)

### Fly dataset information

mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))

tx3fly <- getBM(attributes = c("external_gene_name","ensembl_gene_id","flybase_gene_id","ensembl_transcript_id","transcript_length"),mart=mart)
#len2fly <-getBM(attributes = c("external_gene_name","3_utr_start","3_utr_end"),mart=mart)
#colnames(len2fly) <- c("external_gene_name","UTR_start","UTR_stop")
#len2fly = len2fly %>% replace_na(list(UTR_start=0 ,UTR_stop = 0)) %>% dplyr::mutate(canonical_length=abs(UTR_start-UTR_stop)) %>% 
#  dplyr::group_by(external_gene_name) %>% dplyr::summarise(transcript = n(),UTR_lengths=paste(canonical_length,collapse = ","),canonical_length=max(canonical_length))
#len2fly <- unique(merge(len2fly,tx3fly[,1:2],by="external_gene_name",all.x=T))

#import dataset in GRanges
List_dataset=list()
gtf_file <- "C:/Users/moroa/Desktop/LAB/Genomes/Drosophila/flybase_dmel_r6.21/gtf/dmel-all-r6.21_GB4.gtf"

###Complete Drosophila datasets
Tx_DB <- makeTxDbFromGFF("C:/Users/moroa/Desktop/LAB/Genomes/Drosophila/flybase_dmel_r6.21/gtf/dmel-all-r6.21.gtf" , format = "gtf")
seqlevels(Tx_DB)[which(seqlevels(Tx_DB)=="mitochondrion_genome")]<-"chrM"
seqlevelsStyle(Tx_DB) <- "UCSC"
Tx_DB <- keepSeqlevels(Tx_DB, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY','chrM'))

### Gene dataset
Genes <- GenomicFeatures::genes(Tx_DB)


### Stop dataset
stop_codon<- gtf2gr(gtf_file,gtf_feature=c("stop_codon"),transcript_id="transcript_id",gene_id="gene_id")
seqlevels(stop_codon)[which(seqlevels(stop_codon)=="chrmitochondrion_genome")]<-"chrM"
stop_codon <- keepSeqlevels(stop_codon, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY','chrM'), pruning.mode="coarse")
MAX_stop_codon <- as(data.frame(stop_codon) %>% dplyr::group_by(gene_id) %>% dplyr::summarise(seqnames=dplyr::first(seqnames),strand=dplyr::first(strand), 
                                                                         start=ifelse(strand=="+",max(start,na.rm = F),min(start,na.rm = F)),end=ifelse(strand=="+",max(end),min(end)), 
                                                                         transcript_id=paste0(transcript_id,collapse = ","), exon_id=paste0(exon_id,collapse = ",")) %>% 
  dplyr::mutate(start_UTR=ifelse(strand=="+",end,start)),"GRanges")

MIN_stop_codon <- as(data.frame(stop_codon) %>% dplyr::group_by(gene_id) %>% dplyr::summarise(seqnames=dplyr::first(seqnames),strand=dplyr::first(strand), 
                                                                                              start=ifelse(strand=="+",min(start,na.rm = F),max(start,na.rm = F)),end=ifelse(strand=="+",min(end),max(end)), 
                                                                                              transcript_id=paste0(transcript_id,collapse = ","), exon_id=paste0(exon_id,collapse = ",")) %>% 
                       dplyr::mutate(start_UTR=ifelse(strand=="+",end,start)),"GRanges")

### 3'UTR dataset
threeUTR <- gtf2gr(gtf_file,gtf_feature=c("3UTR"),transcript_id="transcript_id",gene_id="gene_id")
threeUTR <- keepSeqlevels(threeUTR, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY'), pruning.mode="coarse")
threeUTR <- as(merge(data.frame(threeUTR) %>% dplyr::group_by(gene_id) %>% 
                       dplyr::summarise(seqnames=dplyr::first(seqnames),strand=dplyr::first(strand), 
                                        start=ifelse(strand=="+",NA,min(start,na.rm = F)),end=ifelse(strand=="+",max(end),NA), 
                                        transcript_id=paste0(transcript_id,collapse = ","), exon_id=paste0(exon_id,collapse = ",")),data.frame(MAX_stop_codon)[,c(6,9)],by="gene_id") %>% 
        dplyr::mutate(start=ifelse(is.na(start),start_UTR,start),end=ifelse(is.na(end),start_UTR,end),stop_UTR=ifelse(strand=="+",end,start))%>%
        dplyr::mutate(end=(ifelse(end-start<0,start,end))) %>% 
        dplyr::mutate(start=ifelse(strand=="=",start,start-25),end=ifelse(strand=="+",end+25,end)),"GRanges")

### Ambiguos
ambiguous <- gtf2gr(gtf_file,gtf_feature=c("ambiguous"),transcript_id="transcript_id",gene_id="gene_id")
ambiguous <- keepSeqlevels(ambiguous, c('chr2L','chr2R','chr3L','chr3R','chrX'), pruning.mode="coarse")
ambiguous@ranges@NAMES <- as.character(seq(1,length(ambiguous)))
ambiguous <- reduce(ambiguous)

### Extended UTR + 5 kb
Extend_UTR <- as(data.frame(threeUTR) %>% 
                   dplyr::mutate(start=ifelse(strand=="+",stop_UTR,stop_UTR-5000),
                                 end=ifelse(strand=="+",stop_UTR+5000,stop_UTR)) %>%
                   dplyr::select(seqnames,start,end,strand,gene_id,stop_UTR), "GRanges")

### 5'UTR dataset
fiveUTR <- gtf2gr(gtf_file,gtf_feature=c("5UTR"),transcript_id="transcript_id",gene_id="gene_id")
fiveUTR <- keepSeqlevels(fiveUTR, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY'), pruning.mode="coarse")
fiveUTR@ranges@NAMES <- as.character(seq(1,length(fiveUTR)))
fiveUTR <- as(as.data.frame(fiveUTR) %>% dplyr::group_by(gene_id) %>% 
                dplyr::summarise(seqnames=dplyr::first(as.character(seqnames)),start=min(start), end=max(end),
                                 strand=dplyr::first(strand),transcript_n=n(), transcript_id=paste(transcript_id,collapse = ',')) %>% 
                dplyr::arrange(gene_id),"GRanges")

### Exon dataset
Exon <- gtf2gr(gtf_file,gtf_feature=c("exon"),transcript_id="transcript_id",gene_id="gene_id")
seqlevels(Exon)[which(seqlevels(Exon)=="chrmitochondrion_genome")]<-"chrM"
Exon <- keepSeqlevels(Exon, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY','chrM'), pruning.mode="coarse")
Exon@ranges@NAMES <- as.character(seq(1,length(Exon)))
Exon <- as(as.data.frame(Exon) %>% dplyr::group_by(gene_id) %>% 
             dplyr::summarise(seqnames=dplyr::first(as.character(seqnames)),start=min(start), end=max(end),
                              strand=dplyr::first(strand),transcript_n=n(), transcript_id=paste(transcript_id,collapse = ',')) %>% 
             dplyr::arrange(gene_id),"GRanges")

### Intron dataset
introns = unlist(intronsByTranscript(Tx_DB,use.names = T))
mcols(introns)$transcript_id <- names(introns)
names(introns) <- seq(1,152647)
introns <- as(merge(data.frame(introns),tx3fly[,3:4],by.x="transcript_id",by.y="ensembl_transcript_id",all.x=T) %>% dplyr::select(seqnames,start,end,gene_id=flybase_gene_id,transcript_id),"GRanges")

UTR.with.introns <-findOverlaps(threeUTR,introns, type = "any")
list_UTR_with_introns <- unique(data.frame(threeUTR[queryHits(UTR.with.introns)], introns[subjectHits(UTR.with.introns)],stringsAsFactors=F)%>%
                                  separate_rows(transcript_id,sep=",") %>% dplyr::select(gene_id))




### in list
List_dataset["Exon"] = Exon 
List_dataset["Extend_5kb"] = Extend_UTR
List_dataset["fiveUTR"] = fiveUTR
List_dataset["Introns"] = introns
List_dataset["Stop_MAX"] = MAX_stop_codon 
List_dataset["Stop_MIN"] = MIN_stop_codon 
List_dataset["STOP_codon"] = MAX_stop_codon 
List_dataset["threeUTR"] = threeUTR
List_dataset["UTR_with_introns"] = list_UTR_with_introns
List_dataset["Genes"] = Genes 
List_dataset["Ambiguous"] = ambiguous

rm(list=c("Exon","ambiguous","Extend_UTR","fiveUTR","gtf_file","introns","Genes","mart","MAX_stop_codon","MIN_stop_codon","stop_codon","threeUTR","UTR.with.introns","Tx_DB","list_UTR_with_introns"))
