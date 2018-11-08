#load and clean the gtf file obtain by Java script and modify using bash to have chromosome name as chrX
import_atlas <- function(atlas_file){
  atlas <-  rtracklayer::import.gff2(atlas_file)																					#import gtf file
  name = strsplit(strsplit(atlas_file,"/")[[1]][8],"[.]")[[1]][1]
  mcols(atlas) <-  elementMetadata(atlas)[!names(elementMetadata(atlas)) %in% c("cluster", "source", "type", "score", "phase",name)]	#get rid of columns that are not useful
  #mcols(atlas) <- elementMetadata(atlas)[, !(grepl("^n", colnames(elementMetadata(atlas))))]
  mcols(atlas) <- elementMetadata(atlas)[, !(grepl("\\.3p", colnames(elementMetadata(atlas))))]									#eliminate 3p and untrimmed counts
  mcols(atlas) <- elementMetadata(atlas)[, !(grepl("\\.un", colnames(elementMetadata(atlas))))]
  temp_df <- as.data.frame(elementMetadata(atlas))																				# get rid of NAs and substitute with 0
  temp_df <- apply(temp_df, 2, function(x) as.numeric(x))
  temp_df[is.na(temp_df)] <- 0
  colnames(temp_df) <- sub("^X\\d+_(.+)", "\\1", colnames(temp_df))																#change names of the cell lines, get rid of the X number part
  elementMetadata(atlas) <- temp_df
  names(atlas) <- atlas$site_id																									#name cluster by site_id unique identifier
  seqlevels(atlas)[which(seqlevels(atlas)=="mitochondrion_genome")]<-"chrM"
  seqlevelsStyle(atlas) <- "UCSC"
  sl = seqlevels(atlas)[grep("chr", seqlevels(atlas))]
  atlas <- keepSeqlevels(atlas, sl,pruning.mode ="coarse")
  return(atlas)
}

#Without: took one list and deleate the element in common with a second list
"%w/o%" <- function(x, y) x[!x %in% y]

#### fasta of atlas_GR ranges
get_flankseq <- function(atlas_GR, extencion=40 ){
  require(BSgenome.Dmelanogaster.UCSC.dm6)
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  a <- atlas_GR
  ranges(a) <- ranges(a) + extencion
  end_sequence <- getSeq(genome,a)
  return(end_sequence)
}


####
#CountFeature <- function(atlas,gtf="~/Drosophila/dmel_r6.21_FB2018_02/gtf/dmel-all-r6.21.gtf",extend=1000){
#  require(tidyverse)
#  require(GenomicRanges)
#  require(GenomicFeatures)
#  Tx_DB <- makeTxDbFromGFF(gtf , format = "gtf")
#  seqlevelsStyle(Tx_DB) <- "UCSC"
#  Tx_DB <- keepSeqlevels(Tx_DB, c('chr2L','chr2R','chr3L','chr3R','chr4','chrX','chrY'))
#  tot <- length(atlas)
#  sub1 <- atlas[!(atlas %over% threeUTRsByTranscript(Tx_DB))]
#  utr3 <- tot - length(sub1)
#  Extend_1k <- as(as.data.frame(fiveUTRsByTranscript(Tx_DB)) %>% dplyr::mutate(n_start=ifelse(strand=="+",end+1,start-extend-1),n_end=ifelse(strand=="+",end+extend+1,start-1),stop_codon=ifelse(strand=="+",start-1,end+1)) %>% dplyr::select(seqnames,start=n_start,end=n_end,strand,UTR_s=start,UTR_e=end,stop_codon), "GRanges")
#  sub2 <- sub1[!(sub1 %over% Extend_1k)]
#  utr3_ext <- length(sub1) - length(sub2)
#  sub3 <- sub2[!(sub2 %over% fiveUTRsByTranscript(Tx_DB))]
#  utr5 <- length(sub2) - length(sub3)
#  sub4 <- sub3[!(sub3 %over% exonsBy(Tx_DB))]
#  exon <- length(sub3) - length(sub4)
#  sub5 <- sub4[!(sub4 %over% intronsByTranscript(Tx_DB))]
#  intron <- length(sub4) - length(sub5)
#  intragenic <- tot - utr3 - utr3_ext - utr5 - exon - intron
#  dfeature <- data.frame(utr3,utr3_ext,utr5,exon,intron,intragenic)
#  return(dfeature) 
#}
#####
CountFeature_2 <- function(atlas,expand=0,expand_ann=0){
  atlas$type <- as.character("Intragenic")
  Over.exp.1kb <- findOverlaps(atlas,List_dataset$Extend_1kb, type = "within")
  atlas[(unique(queryHits(Over.exp.1kb)))]$type <- "CDS+1kb"
  Over.exon <- findOverlaps(atlas,List_dataset$Exon, type = "within")
  atlas[(unique(queryHits(Over.exon)))]$type <- "Exons"
  Over.5UTRs <- findOverlaps(atlas,List_dataset$fiveUTR, type = "within")
  atlas[(unique(queryHits(Over.5UTRs)))]$type <- "5UTRs"
  Over.CDS.rdc <- findOverlaps(atlas,List_dataset$CDS, type = "within")
  atlas[(unique(queryHits(Over.CDS.rdc)))]$type <- "CDSs"
  Over.introns <- findOverlaps(atlas,List_dataset$Introns, type = "within")
  atlas[(unique(queryHits(Over.introns)))]$type <- "Introns"
  Over.exp.Piero <- findOverlaps(atlas,List_dataset$threeUTR_Piero+expand, type = "within")
  atlas[(unique(queryHits(Over.exp.Piero)))]$type <- "3UTRs_confidence"
  Over.all.3UTRs <- findOverlaps(atlas,List_dataset$threeUTR+expand_ann, type = "within")
  atlas[(unique(queryHits(Over.all.3UTRs)))]$type <- "3UTRs"
  atlas[unique(queryHits(Over.exon[queryHits(Over.exon) %>% duplicated]))]$type <- "ambiguous"
  return(atlas)
}
#####
CountFeature_3 <- function(atlas,expand=0,expand_ann=0){
  atlas$type <- as.character("Intragenic")
  atlas$start_UTR <- NA
  atlas$stop_UTR <- NA
  atlas$gene_id <- NA
  Over.exp.5kb <- findOverlaps(atlas,List_dataset$Extend_5kb, type = "any")
  atlas[queryHits(Over.exp.5kb)]$type <- "5kb"
  atlas[queryHits(Over.exp.5kb)]$stop_UTR <- List_dataset$Extend_5kb[subjectHits(Over.exp.5kb)]$stop_UTR
  atlas[queryHits(Over.exp.5kb)]$gene_id <- List_dataset$Extend_5kb[subjectHits(Over.exp.5kb)]$gene_id
  Over.exon <- findOverlaps(atlas,List_dataset$Exon, type = "any")
  atlas[queryHits(Over.exon)]$type <- "Exons"
  atlas[queryHits(Over.exon)]$gene_id <- List_dataset$Exon[subjectHits(Over.exon)]$gene_id
  Over.introns <- findOverlaps(atlas,List_dataset$Introns, type = "any")
  atlas[queryHits(Over.introns)]$type <- "Introns"
  atlas[queryHits(Over.introns)]$gene_id <- List_dataset$Introns[subjectHits(Over.introns)]$gene_id
  Over.5UTRs <- findOverlaps(atlas,List_dataset$fiveUTR, type = "any")
  atlas[queryHits(Over.5UTRs)]$type <- "5UTRs"
  atlas[queryHits(Over.5UTRs)]$gene_id <- List_dataset$fiveUTR[subjectHits(Over.5UTRs)]$gene_id
  Over.all.3UTRs <- findOverlaps(atlas,List_dataset$threeUTR, type = "any")
  atlas[queryHits(Over.all.3UTRs)]$type <- "3UTRs"
  atlas[queryHits(Over.all.3UTRs)]$gene_id <- List_dataset$threeUTR[subjectHits(Over.all.3UTRs)]$gene_id
  atlas[queryHits(Over.all.3UTRs)]$start_UTR <- List_dataset$threeUTR[subjectHits(Over.all.3UTRs)]$start_UTR 
  atlas[queryHits(Over.all.3UTRs)]$stop_UTR <- List_dataset$threeUTR[subjectHits(Over.all.3UTRs)]$stop_UTR
  Over.ambiguous <- findOverlaps(atlas,List_dataset$Ambiguous, type = "any")
  atlas[queryHits(Over.ambiguous)]$type <- "ambiguous"
  atlas[unique(queryHits(Over.all.3UTRs[queryHits(Over.all.3UTRs) %>% duplicated]))]$type <- "ambiguous2"
  return(atlas)
}

#### Find # introns in UTR
find_intron_UTR <- function(x){
  x$intron_UTR <- 0
  number <- x$exon_rank %>% length()
  if (number > 1){x$intron_UTR <- number}
  return(x)
}

#### merge UTR with multiple introns
merge_intron_UTR <- function(x){
  if ((x$intron_UTR %>% length) > 1){
    intron_UTR <- x$intron_UTR[1]
    tx_name <- x$tx_name[1]
    gene_id <- x$gene_id[1]
    x <- GRanges(seqnames = seqnames(x)[1],ranges =  IRanges(start = min(start(x)), end = max(end(x))), strand = as.vector(strand(x))[1])
    x$tx_name <- tx_name
    x$gene_id <- gene_id
    x$intron_UTR <- intron_UTR
    return(x)
  } else{return(x)}
}


add_PAS_GR <- function(GR, PAS_list, genome){
  length_GR_col <- (GR %>% mcols %>% length)+1
  GR <- assign_major_PAS(PAS_set = major_UTR_PAS, GR = GR, genome = genome)
  print(GR)
  GR$PAS <- "no_major_PAS"
  for (i in 10:1){
    a <- as.character(PAS_list[i])
    GR[mcols(GR)[,a] > 0]$PAS <- a
  }
  GR$PAS_number <- apply(as.data.frame(mcols(GR))[length_GR_col:(length_GR_col+9)], 1, sum)
  GR$PAS_labels <- GR$PAS
  GR[GR$PAS_labels != "AATAAA" & GR$PAS_labels != "ATTAAA" & GR$PAS_labels != "AATATA" & GR$PAS_labels != "no_major_PAS"]$PAS_labels <- "other_PAS"
  GR_df <- as.data.frame(mcols(GR))
  return(GR_df)
}


my_matrix_new <- function(atlas_ends, genome){
  a <- atlas_ends
  ranges(a) <- ranges(a) + 100L
  a <- getSeq(genome,a)
  a <- consensusMatrix(a, as.prob = TRUE, baseOnly=TRUE)
  a <- a[1:4,]
  a <- as.data.frame(a)
  names(a) <- -100:100
  a <- as.data.frame(a)
  a$nucleotide <- rownames(a)
  a <- tidyr::gather(a, position, percentage, 1:200)
  ggplot(a, aes(x=as.numeric(position), y=percentage, colour = nucleotide)) +
    geom_line(size=1.1) +
    theme_bw() +
    theme(legend.key = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.title=element_blank(),
          legend.title =element_blank(),
          legend.position = c(.92, .82),
          legend.text = element_text(size = 20),
          plot.margin = unit(c(1,1,1,1), "cm"),
          legend.key.size = unit(0.8, "cm")) +
    scale_x_continuous(expand = c(0, 0),
                       breaks=seq(-100,100,25), limits= c(-100,100)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0,1)) +
    scale_color_manual(labels = c("A", "C", "G", "U"), values=c("red2", "goldenrod2", "limegreen", "blue3"))
}
