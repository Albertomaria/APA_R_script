setwd("C:/Users/moroa/Desktop/3seq_analysis")

## This script use a gtf to associate gene, genomic feature, Stop_codon, termination of 3'UTR to atlas generate with Sol Script


#Libraries
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(TransView)
library(IRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(scales)
library(edgeR)
library(matrixStats)

########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
### load function, dataset and save palet for graphs                                                      ###
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

source("C:/Users/moroa/Desktop/R_Script/Functions.R")
source("C:/Users/moroa/Desktop/R_Script/Dataset_from_GB.R")
palet_feature <- c("#9f2f42","#66c2a5","#5e4fa2", "#3288bd", "#ef6c42","#f6ad60","#fde08b","#bebebe")  #d5404f -> red
palet_euler <- c("#d62d20", "#008744","#0057e7","#ffa700","#a700ff","#00ffa7","#bebebe")
palet_boxplot <- c("#BEBEBE","#85b7d7","#296c97","#e58a95","#bf3647","#f9f9ea","#d8d866","#bebe00","#a5d09a","#72aa63")


########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
### load atlas and combine atlas in one master atlas                                                      ###
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

gtf <- list.files("E:/3Seq/atlas_gtf/")

all_atlas <- list()
for (file in gtf[grep(".gtf",gtf)]){
  print(file)
  name=gsub("[[:punct:]]","-",strsplit(file,"\\.")[[1]][1])
  gtf_df=import_atlas(atlas_file = paste0("E:/3Seq/atlas_gtf/",file))
  all_atlas[name]<-gtf_df
}
rm(list = c("name","gtf_df"))

who <- names(all_atlas[grep("S2|Head|L1|Lijuan",names(all_atlas))])
who <- who[!grepl("MutR1|Nuc|Cyt|CME",who)]

atlas <- all_atlas[who]

atlas_df <- data.frame(termination=NA)
cnames <- "termination"

for (n in who){
  df <- data.frame(atlas[[n]])
  df <- df %>% dplyr::mutate(termination = paste(seqnames,start,strand, sep = ":")) %>% dplyr::select(termination,n)
  atlas_df <- merge(atlas_df,df,by="termination",all=T)
  cnames <- c(cnames,n)
}
colnames(atlas_df) <- cnames
rm(list=c("cnames","who","df"))

atlas_df[is.na(atlas_df)] <- 0
atlas_total <- as(atlas_df %>% dplyr::filter(termination!=0) %>% tidyr::separate(termination,into=c("seqnames","start","strand"),sep=":",remove=F) %>% dplyr::mutate(end=start),"GRanges")


########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
### rename colum and caluclate CPM                                                                        ###
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

atlas_df <- data.frame(atlas_total) %>% dplyr::select(seqnames,start,end,strand,termination,dplyr::everything())
colnames(atlas_df) <- c("seqnames","start","end","strand","termination","width",
                        "Control_S2.1","Control_S2.2","Control_S2.3",
                        "L1_Elav_KO.1","L1_Elav_KO.2","L1_Elav_KO.3","L1_Elav_Rescue.1","L1_Elav_Rescue.2","L1_Elav_Rescue.3",
                        "Elav_Mut3X.1","Elav_Mut3X.2","Elav_Mut3X.3","Elav_WT.1","Elav_WT.2","Elav_WT.3",
                        "Fne_Mut3X.1","Fne_WT.1","Fne_WT.2","Fne_WT.3","Head.1",
                        "ime4_Dcat","ime4_Plus","mettl14_Df","mettl14_Plus", 
                        "Rbp9_Mut3X.1","Rbp9_WT.1","Rbp9_WT.2","Rbp9_WT.3")
atlas_RPM <- cbind(atlas_df[,1:6],cpm(atlas_df[,7:ncol(atlas_df)]))
atlas_RPM_GR <- as(atlas_RPM,"GRanges")

rm(list=c("file","gtf","n"))


########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
### grouping the termination in cluster of 25 nt aroung the termination (50nt in total)                   ###
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

grouping <- as(data.frame(reduce(atlas_RPM_GR+25)) %>% dplyr::mutate(termination = paste0(seqnames,":",start,"-",end,":",strand)),"GRanges")
new_group <- findOverlaps(atlas_RPM_GR,grouping)

temp <- atlas_RPM_GR

temp$cluster <- NA
temp[(queryHits(new_group))]$cluster <- grouping[(subjectHits(new_group))]$termination

atlas_RPM_group <- merge(data.frame(temp) %>% dplyr::group_by(cluster) %>% dplyr::summarise_at(vars(dplyr::matches("Control|Head|Elav|Fne|Rbp9|ime4|mettl14")),funs(max)),
                         data.frame(temp) %>% dplyr::group_by(cluster) %>% dplyr::summarise(termination=paste0(termination,collapse = ","),
                                                                                                           n_termination=n()),by="cluster") %>%
  tidyr::separate(cluster,into=c("seqnames","ranges","strand"),sep=":",remove=F) %>% tidyr::separate(ranges,into=c("start","end"),sep="-",convert=T) %>%
  dplyr::select(seqnames,start,end,strand,cluster,termination,n_termination, everything())

atlas_RPM_group_GR <- as(atlas_RPM_group,"GRanges")

#atlas_RPM_group <- data.frame(atlas_RPM_group_GR) %>% dplyr::select(seqnames,start,end,width,strand,gene_id,type,start_UTR,stop_UTR,cluster,termination,n_termination, everything())

rm(list = c("new_group","grouping"))

#atlas_RPM_filter <- cbind(atlas_RPM_group,data.frame("MAX_RPM"=rowMaxs(as.matrix(atlas_RPM_group[,13:ncol(atlas_RPM_group)])))) %>% dplyr::filter(MAX_RPM>5)
#atlas_RPM_filter_GR <- as(atlas_RPM_filter,"GRanges")
#atlas_RPM_filter_GR <- CountFeature_3(atlas_RPM_filter_GR)
#atlas_RPM <- data.frame(CountFeature_3(atlas_RPM_GR))
atlas_RPM_group_GR <- data.frame(CountFeature_3(atlas_RPM_group_GR))
#atlas_RPM_filter <- data.frame(atlas_RPM_filter_GR)

table(data.frame(atlas_RPM_group_GR)$type)

Distribution_5kb_group <- data.frame(atlas_RPM_group_GR) %>% dplyr::mutate (distance=ifelse(strand=="+" & type=="5kb",end-stop_UTR,ifelse(strand=="-" & type=="5kb",stop_UTR-start,NA)),
                                                                                            Breaks = cut(distance, breaks = c(seq(0, 50, 25),100,500,1000,5000),na.rm=T))

  
write.table(atlas_RPM_filter,"termination_distribution.tab",sep="\t",quote = F, row.names = F)


write.table(Distribution_5kb_group,"termination_distribution_5kb_Count.tab",sep="\t",quote = F, row.names = F)

##### graphical distribution
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
Feature_distribution <- data.frame(table(atlas_RPM$type))  
Feature_distribution_filtered <- data.frame(table(atlas_RPM_filter$type))
Feature_distribution_grouped <- data.frame(table(atlas_RPM_group$type)) 
Feature_distribution <- cbind(Feature_distribution,Feature_distribution_filtered[,2],Feature_distribution_grouped[,2])
colnames(Feature_distribution) <- c("Feature",paste0("Total: # ",dim(atlas_RPM)[[1]]),paste0("Filter_mean.5: # ",dim(atlas_RPM_filter)[[1]]),paste0("Grouped_10nt: # ",dim(atlas_RPM_group)[[1]]))
Feature_distribution <- Feature_distribution %>% tidyr::gather(key="Condition",value="Reads",-Feature) %>% dplyr::group_by(Condition) %>% dplyr::mutate(Percentage=Reads/sum(Reads)*100)


Feature_distribution$Feature <- factor(Feature_distribution$Feature,levels = c("3UTRs","5UTRs","Exons","Introns","5kb","Intragenic","ambiguous"))
Feature_distribution$Condition <- factor(Feature_distribution$Condition,levels = c(paste0("Total: # ",dim(atlas_RPM)[[1]]),paste0("Filter_mean.5: # ",dim(atlas_RPM_filter)[[1]]),
                                                                                   paste0("Grouped_10nt: # ",dim(atlas_RPM_group)[[1]])))
tiff("New_grap_per_presentation/Distribution_feature_3.tif",width = 1000, height = 1000)
ggplot(Feature_distribution, aes(x=Feature,y=Reads,fill=Feature)) + geom_bar(stat="identity") +
  scale_fill_manual(values = palet_feature) + labs(x="", y="Reads", fill="Genomic Feature") +
  facet_wrap(~Condition) + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(size=24,angle=45,hjust = 1),
        axis.text.y = element_text(size=24), axis.title.y = element_text(size=28,face="bold"),
        strip.text.x = element_text(size=16,face="bold"),
        strip.background = element_rect(fill="#dcdcc6"),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20))
dev.off()
rm(list=c("Feature_distribution","Feature_distribution_filtered","Feature_distribution_grouped"))

##### graphical distribution 5kb
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
Distribution_5kb_total <- atlas_RPM %>% dplyr::filter(type=="5kb") %>% dplyr::mutate (distance=ifelse(strand=="+",end-stop_UTR,stop_UTR-start),
                                                                                      Bricks = cut(distance, breaks = c(seq(0, 50, 25),100,500,1000,5000)),
                                                                                      RPM_mean = rowMeans(.[,10:23])) %>%
  dplyr::group_by(Bricks) %>% dplyr::summarise(Reads=n(),RPM_mean=mean(RPM_mean)) %>% dplyr::mutate (Sample=paste0("Total: # ",sum(Reads)))

Distribution_5kb_filter <- atlas_RPM_filter %>% dplyr::filter(type=="5kb") %>% dplyr::mutate (distance=ifelse(strand=="+",end-stop_UTR,stop_UTR-start),
                                                                                              Bricks = cut(distance, breaks = c(seq(0, 50, 25),100,500,1000,5000))) %>%
  dplyr::group_by(Bricks) %>% dplyr::summarise(Reads=n(),RPM_mean=mean(MAX_RPM)) %>% dplyr::mutate (Sample=paste0("Filter_mean.5: # ",sum(Reads)))

Distribution_5kb_group <- atlas_RPM_group %>% dplyr::filter(type=="5kb") %>% dplyr::mutate (distance=ifelse(strand=="+",end-stop_UTR,stop_UTR-start),
                                                                                            Bricks = cut(distance, breaks = c(seq(0, 50, 25),100,500,1000,5000)),
                                                                                            RPM_mean = rowMeans(.[,9:36])) %>%
  dplyr::group_by(Bricks) %>% dplyr::summarise(Reads=n(),RPM_mean=mean(RPM_mean)) %>% dplyr::mutate (Sample=paste0("Grouped_10nt: # ",sum(Reads))) %>% dplyr::filter(!is.na(Bricks))

Distribution_5kb <- rbind(Distribution_5kb_total,Distribution_5kb_filter,Distribution_5kb_group) %>% dplyr::group_by(Sample) %>% dplyr::mutate(Percentage=Reads/sum(Reads))
Distribution_5kb$Sample <- factor(Distribution_5kb$Sample,levels = c("Total: # 2063","Filter_mean.5: # 368","Grouped_10nt: # 1539"))

tiff("New_grap_per_presentation/Distribution_5k_3.tif",width = 1000, height = 1000)
ggplot(Distribution_5kb, aes(x=Bricks,y=Reads,fill=RPM_mean)) + geom_bar(stat="identity",color="black") +
  labs(x="", y="Reads", fill="average(RPM)") + scale_fill_gradient(low = "#FFFBEF", high = "#ef6c42") +
  facet_wrap(.~Sample) + 
  scale_x_discrete(labels=c("(0,25]" = "0-25", "(25,50]"="26-50","(50,100]"="51-100","(100,500]"="101-500","(500,1e+03]"="501-1000", "(1e+03,5e+03]"="+ 1000")) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(size=24,angle=45,hjust = 1),
        axis.text.y = element_text(size=24), axis.title.y = element_text(size=28,face="bold"),
        strip.text.x = element_text(size=16,face="bold"),
        strip.background = element_rect(fill="#dcdcc6"),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20))
???dev.off()

rm(list=c("Distribution_5kb","Distribution_5kb_total","Distribution_5kb_filter","Distribution_5kb_group"))
