library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
library(readr)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(TransView)

dm6_alias <- read_delim("C:/Users/moroa/igv/genomes/dm6_alias.tab", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


colnames(dm6_alias) <- c("short","long")
gtf_file <- "C:/Users/moroa/Desktop/LAB/Genomes/Drosophila/flybase_dmel_r6.21/gtf/dmel-all-r6.21.gtf"
gtf_Piero <- rtracklayer::import('C:/Users/moroa/Desktop/LAB/GB_Drosophila_landscape/total_Drosophila_from4.gtf')

stop_codon<- gtf2gr(gtf_file,gtf_feature=c("stop_codon"),transcript_id="transcript_id",gene_id="gene_id")
seqlevels(stop_codon)[which(seqlevels(stop_codon)=="chrmitochondrion_genome")]<-"chrM"
stop_codon <- keepSeqlevels(stop_codon, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY','chrM'), pruning.mode="coarse")
MAX_stop_codon <- as(data.frame(stop_codon) %>% dplyr::group_by(gene_id) %>% dplyr::summarise(seqnames=dplyr::first(seqnames),strand=dplyr::first(strand), 
                                                                                              start=ifelse(strand=="+",max(start,na.rm = F),min(start,na.rm = F)),end=ifelse(strand=="+",max(end),min(end)), 
                                                                                              transcript_id=paste0(transcript_id,collapse = ","), exon_id=paste0(exon_id,collapse = ",")) %>% 
                       dplyr::mutate(fine=ifelse(strand=="+",end,start)),"GRanges")

knowledge <- data.frame("unknown" = c("FBgn0011241","FBgn0015541","FBgn0016700","FBgn0019886","FBgn0019947","FBgn0021818","FBgn0024332","FBgn0025286","FBgn0025687","FBgn0028665","FBgn0029128",
                                      "FBgn0029173","FBgn0031781","FBgn0032679","FBgn0032915","FBgn0033059","FBgn0033342","FBgn0034133","FBgn0035707","FBgn0037645","FBgn0038035","FBgn0038220",
                                      "FBgn0038341","FBgn0039177","FBgn0042650","FBgn0050007","FBgn0050118","FBgn0051728","FBgn0053123","FBgn0053139","FBgn0086901","FBgn0250876","FBgn0260441",
                                      "FBgn0260460","FBgn0261241","FBgn0261548","FBgn0261593","FBgn0261613","FBgn0262735","FBgn0262742","FBgn0267911","FBgn0000497","FBgn0000054","FBgn0004009",
                                      "FBgn0000556","FBgn0005592","FBgn0001974","FBgn0000409","FBgn0010575","FBgn0003887","FBgn0003732","FBgn0002607"),
                        "known" = c("FBgn0285962","FBgn0285963","FBgn0285937","FBgn0284252","FBgn0284421","FBgn0286070","FBgn0284442","FBgn0285949","FBgn0286198","FBgn0285910","FBgn0284435",
                                    "FBgn0285971","FBgn0284255","FBgn0284256","FBgn0284251","FBgn0284246","FBgn0284436","FBgn0285943","FBgn0286051","FBgn0286204","FBgn0286075","FBgn0285913",
                                    "FBgn0286071","FBgn0286199","FBgn0285879","FBgn0285892","FBgn0284257","FBgn0284244","FBgn0284253","FBgn0284254","FBgn0285955","FBgn0284221","FBgn0286213",
                                    "FBgn0285951","FBgn0285911","FBgn0283741","FBgn0285947","FBgn0284250","FBgn0285926","FBgn0285925","FBgn0284408","FBgn0284247","FBgn0284249","FBgn0284084",
                                    "FBgn0284245","FBgn0285896","FBgn0286027","FBgn0284248","FBgn0285917","FBgn0284243","FBgn0284220","FBgn0285950"))

gtf_Piero_clean <- data.frame(gtf_Piero) %>% dplyr::mutate(type=ifelse(source=="Drosophila_3seq",as.character(genomic_feature),as.character(type)),
                                                           gene_id=ifelse(source=="Drosophila_3seq",FlyBase_ID,gene_id),
                                                           gene_symbol=ifelse(source=="Drosophila_3seq",symbol,gene_symbol),
                                                           Note=ifelse(source=="Drosophila_3seq","major_end",`X.`)) %>%
  dplyr::select(seqnames,start,end,width,strand,source,type,score,phase,gene_id,gene_symbol,transcript_id,transcript_symbol,Note) %>% 
  tidyr::separate_rows(gene_id,sep="_") %>% 
  merge(knowledge,.,by.x="unknown",by.y="gene_id",all.y=T) %>% 
  dplyr::mutate(gene_id=ifelse(!is.na(known),as.character(known),as.character(unknown))) %>% dplyr::select(-known,-unknown) %>%
  merge(dm6_alias,.,by.x="short",by.y="seqnames",all.y=T) %>%
  dplyr::mutate(seqnames=ifelse(is.na(long),as.character(short),as.character(long))) %>% dplyr::select(-short,-long) %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(gene_symbol=dplyr::first(na.omit(gene_symbol))) %>% 
  dplyr::mutate(Note=ifelse(is.na(Note),".",Note))

gtf_Piero_clean_1 <-  merge(data.frame(MAX_stop_codon)[,c(6,9)],gtf_Piero_clean, by="gene_id",all.y=T) %>% 
  dplyr::mutate(start=ifelse(Note=="major_end" & strand=="+" & !is.na(fine),fine,start)) %>% 
  dplyr::mutate(end=ifelse(Note=="major_end"&strand=="-" & !is.na(fine),fine,end)) %>% 
  dplyr::mutate(Note=ifelse(Note=="major_end",paste0("ME_ext_",end-start),Note)) %>%
  dplyr::mutate(newV1=ifelse(start>end,end,start),newV2=ifelse(start>end,start,end),start=newV1,end=newV2) %>% dplyr::select(-fine,-newV1, -newV2) 

gtf_Piero_clean_1[is.na(gtf_Piero_clean_1)] <- ""

View(gtf_Piero_clean_1)

gtf_Piero_GR <- as(gtf_Piero_clean_1,"GRanges")

rtracklayer::export.gff2(gtf_Piero_clean_1,"C:/Users/moroa/Desktop/LAB/Genomes/Drosophila/flybase_dmel_r6.21/gtf/dmel-all-r6.21_GB4.gtf")
