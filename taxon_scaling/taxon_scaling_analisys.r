# Script used for the normalisation of the transcriptomic data. It is adapted from  Klingenberg H., Meinicke P. 2017 and Christel S. et al. 2018


library(tidyverse)
library(DESeq2)
library(stringr)
library(biobroom)
library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)
library(grid)

#Load data
#Load sample information tableLoad featurecounts table

#set working directory relative to repository
setwd=here::here()
samples_Mono_Meta <- read_tsv("featureCounts/samples.txt")

#load feature count files
#bacteroides_a
input_files_BACA <- list.files("featureCounts/", 
                          pattern='bacteroides_a', full.names=TRUE)
input_files_BACA %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_BACA

#bacteroides_b
input_files_BACB <- list.files("featureCounts/", 
                          pattern='bacteroides_b', full.names=TRUE)
input_files_BACB %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_BACB

#bacteroides_c
input_files_BACC <- list.files("featureCounts/", 
                          pattern='bacteroides_c', full.names=TRUE)
input_files_BACC %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_BACC

#citrobacter
input_files_CITR <- list.files("featureCounts/", 
                          pattern='citrobacter', full.names=TRUE)
input_files_CITR %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_CITR

#fusobacterium
input_files_FUS <- list.files("featureCounts/", 
                          pattern='fusobacterium', full.names=TRUE)
input_files_FUS %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_FUS

#kerstersia
input_files_KER <- list.files("featureCounts/", 
                          pattern='kerstersia', full.names=TRUE)
input_files_KER %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_KER

#monocercomonoides
input_files_MONO <- list.files("featureCounts/", 
                          pattern='monocercomonoides', full.names=TRUE)
input_files_MONO %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_MONO

#parabacteroides
input_files_PARBA <- list.files("featureCounts/", 
                          pattern='parabacteroides', full.names=TRUE)
input_files_PARBA %>%
  map(read_tsv,comment="#",col_types=cols()) %>%
  set_names(str_extract(.,"RD\\dR\\d")) %>%
  plyr::join_all() -> featcounts_PARBA

#change sample names to meaningful IDs

#bacteroides_a
colnames(featcounts_BACA)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_BACA))
colnames(featcounts_BACA)=plyr::mapvalues(names(featcounts_BACA),
                                     from = samples_Mono_Meta$`Sample ID`,
                                     to = as.character(samples_Mono_Meta$Day))
featcounts_BACA %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_BACA

#bacteroides_b
colnames(featcounts_BACB)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_BACB))
colnames(featcounts_BACB)=plyr::mapvalues(names(featcounts_BACB),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_BACB %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_BACB

#bacteroides_c
colnames(featcounts_BACC)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_BACC))
colnames(featcounts_BACC)=plyr::mapvalues(names(featcounts_BACC),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_BACC %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_BACC

#citrobacter
colnames(featcounts_CITR)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_CITR))
colnames(featcounts_CITR)=plyr::mapvalues(names(featcounts_CITR),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_CITR %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_CITR

#fusobacterium
colnames(featcounts_FUS)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_FUS))
colnames(featcounts_FUS)=plyr::mapvalues(names(featcounts_FUS),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_FUS %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_FUS

#kerstersia
colnames(featcounts_KER)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_KER))
colnames(featcounts_KER)=plyr::mapvalues(names(featcounts_KER),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_KER %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_KER

#monocercomonoides
colnames(featcounts_MONO)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_MONO))
colnames(featcounts_MONO)=plyr::mapvalues(names(featcounts_MONO),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_MONO %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_MONO

#parabacteroides_a
colnames(featcounts_PARBA)=gsub(".+(RD\\dR\\d).+","\\1",colnames(featcounts_PARBA))
colnames(featcounts_PARBA)=plyr::mapvalues(names(featcounts_PARBA),
                                      from = samples_Mono_Meta$`Sample ID`,
                                      to = as.character(samples_Mono_Meta$Day))
featcounts_PARBA %>% gather(sampleID,readcounts,-Geneid,-Chr,-Start,-End,-Strand,-Length) %>%
  mutate(day=gsub("(D\\d)R\\d","\\1",sampleID)) -> featcounts_PARBA

#Normalize feature counts with Deseq2
DESeq2.norm.mat <- function(Xmat,cond,type)
{
  #Xmat = raw count data for one species
  Xmat.col = ncol(Xmat)
  Xmat.row = nrow(Xmat)
  colData <- data.frame(condition = cond, type = type)
  Xmat = round(Xmat)
  storage.mode(Xmat) <- 'integer'
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)
  colData(dds)$condition = factor(colData(dds)$condition, levels=unique(cond))
  dds <- DESeq(dds, quiet = TRUE)
  #normalize the data
  YMat <- Xmat/rep(dds@colData@listData$sizeFactor, each = (Xmat.row))
  return(YMat)
}

DESeq2.result <- function(Xmat,cond,type,CONTRAST)
{
  Xmat.col = ncol(Xmat)
  Xmat.row = nrow(Xmat)
  Xmat = round(Xmat)
  storage.mode(Xmat) <- 'integer'
  colData <- data.frame(condition = cond, type = type)
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)
  
  colData(dds)$condition = factor(colData(dds)$condition,
                                  levels=unique(cond))
  
  #stop DESeq2 from performing additional normalization
  normFactors <- matrix(1,ncol = Xmat.col, nrow = Xmat.row)
  normalizationFactors(dds) <- normFactors
  dds <- DESeq(dds, quiet = F)
  res <- results(dds,contrast=CONTRAST)
  return(res)
}

#conditions based on organism combination
cond.vec=gsub("(D\\d)R\\d","\\1",samples_Mono_Meta$Day)
type.vec=rep("pe",9)

#matrix for BACA
BACA=featcounts_BACA %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(BACA)=BACA$Geneid
BACA=as.matrix(BACA[,-1]  )

#matrix for BACB
BACB=featcounts_BACB %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(BACB)=BACB$Geneid
BACB=as.matrix(BACB[,-1]  )

#matrix for BACC
BACC=featcounts_BACC %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(BACC)=BACC$Geneid
BACC=as.matrix(BACC[,-1]  )

#matrix for CITR
CITR=featcounts_CITR %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(CITR)=CITR$Geneid
CITR=as.matrix(CITR[,-1]  )

#matrix for FUS
FUS=featcounts_FUS %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(FUS)=FUS$Geneid
FUS=as.matrix(FUS[,-1]  )

#matrix for KER
KER=featcounts_KER %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(KER)=KER$Geneid
KER=as.matrix(KER[,-1]  )

#matrix for MONO
MONO=featcounts_MONO %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(MONO)=MONO$Geneid
MONO=as.matrix(MONO[,-1]  )

#matrix for PARBA
PARBA=featcounts_PARBA %>% 
  select(Geneid,sampleID,readcounts) %>%
  spread(sampleID,readcounts) 
rownames(PARBA)=PARBA$Geneid
PARBA=as.matrix(PARBA[,-1]  )

#normalize matrices
BACA_nor=DESeq2.norm.mat(BACA,cond.vec,type.vec)
BACB_nor=DESeq2.norm.mat(BACB,cond.vec,type.vec)
BACC_nor=DESeq2.norm.mat(BACC,cond.vec,type.vec)
CITR_nor=DESeq2.norm.mat(CITR,cond.vec,type.vec)
FUS_nor=DESeq2.norm.mat(FUS,cond.vec,type.vec)
KER_nor=DESeq2.norm.mat(KER,cond.vec,type.vec)
MONO_nor=DESeq2.norm.mat(MONO,cond.vec,type.vec)
PARBA_nor=DESeq2.norm.mat(PARBA,cond.vec,type.vec)

#run deseq for possible combinations

#bacteroides_a
BACA_res1=DESeq2.result(BACA_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
BACA_res2=DESeq2.result(BACA_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
BACA_res3=DESeq2.result(BACA_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#bacteroides_b
BACB_res1=DESeq2.result(BACB_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
BACB_res2=DESeq2.result(BACB_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
BACB_res3=DESeq2.result(BACB_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#bacteroides_c
BACC_res1=DESeq2.result(BACC_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
BACC_res2=DESeq2.result(BACC_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
BACC_res3=DESeq2.result(BACC_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#citrobacter
CITR_res1=DESeq2.result(CITR_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
CITR_res2=DESeq2.result(CITR_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
CITR_res3=DESeq2.result(CITR_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#fusobacterium
FUS_res1=DESeq2.result(FUS_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
FUS_res2=DESeq2.result(FUS_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
FUS_res3=DESeq2.result(FUS_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#kerstersia
KER_res1=DESeq2.result(KER_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
KER_res2=DESeq2.result(KER_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
KER_res3=DESeq2.result(KER_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#monocercomonoides
MONO_res1=DESeq2.result(MONO_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
MONO_res2=DESeq2.result(MONO_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
MONO_res3=DESeq2.result(MONO_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

#parabacteroides_a
PARBA_res1=DESeq2.result(PARBA_nor,cond.vec,type.vec,CONTRAST=c("condition","D3","D2")) %>% tidy %>%mutate(term="D2_D3")
PARBA_res2=DESeq2.result(PARBA_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D3")) %>% tidy %>%mutate(term="D3_D5")
PARBA_res3=DESeq2.result(PARBA_nor,cond.vec,type.vec,CONTRAST=c("condition","D5","D2")) %>% tidy %>%mutate(term="D2_D5")

##create master chief table
master=bind_rows(BACA_res1,BACA_res2,BACA_res3,BACB_res1,BACB_res2,BACB_res3,BACC_res1,BACC_res2,BACC_res3,CITR_res1,CITR_res2,CITR_res3,FUS_res1,FUS_res2,FUS_res3,KER_res1,KER_res2,KER_res3,MONO_res1,MONO_res2,MONO_res3,PARBA_res1,PARBA_res2,PARBA_res3)

##write out table (The identification of DE genes based on FoldChange and p.value is encoded here)
master %>% rename(estimate="log2FoldChange") %>%
  mutate(org=ifelse(grepl("^BACA",gene),"BACA",ifelse(grepl("^BACB",gene),"BACB",ifelse(grepl("^BACC",gene),"BACC",ifelse(grepl("^CITR",gene),"CITR",ifelse(grepl("^FUS",gene),"FUS",ifelse(grepl("^KER",gene),"KER",ifelse(grepl("MONOS",gene),"MONO","PARBA")))))))) %>%
  mutate (diffexpressed=ifelse(log2FoldChange > 1.0 & p.adjusted < 0.05,"UP",ifelse(log2FoldChange < -1.0 & p.adjusted < 0.05,"DOWN", "NO"))) -> master.annot
write.table(master.annot,"featureCounts/comparisons_table.tsv",quote=F,row.names=F,col.names = T,sep="\t")

#Table of DE Genes
#some filters for easier searching
#D3 vs D2
master.selection <- function(sp, day) {
  master.annot %>%
    filter(org==sp) %>% #BACA, BACB, BACC, CITR, FUS, KER, MONO, PARBA
    filter(term==day) %>%
    arrange(desc(abs(log2FoldChange))) -> master.selection.function
  return(master.selection.function)

} 

master.annot %>%
  filter(abs(log2FoldChange)>1) %>%
  filter(p.adjusted < 0.05) %>%
  arrange(desc(log2FoldChange)) -> master.annot.sort

write.table(master.annot.sort,"featureCounts/comparisons_table_sorted.tsv", quote=F,row.names=F,col.names = T,sep="\t")

master.selection(sp = "PARBA", day = "D2_D3") -> PARBA_selection_1
master.selection(sp = "PARBA", day = "D3_D5") -> PARBA_selection_2
master.selection(sp = "PARBA", day = "D2_D5") -> PARBA_selection_3
#BACA, BACB, BACC, CITR, FUS, KER, MONO, PARBA


#Volcano pots
mycolors <- c("royalblue", "grey30", "red2")
names(mycolors) <- c("UP", "NO", "DOWN")

volcano.plot <- function(table, title){
  p <-ggplot(data = table, aes(x=log2FoldChange, y=-log10(p.adjusted), col=diffexpressed)) +
        labs(title = title) +
        geom_point() +
        theme_classic() +
        #geom_text() +
        #geom_text_repel() +
        geom_vline(xintercept =c(-1,1), col="red", linetype = 2 ) +
        geom_hline(yintercept=-log10(0.05),col="red", linetype = 2) +
        scale_color_manual(values= mycolors)
  
  return(p)
}

volcano.plot(table = PARBA_selection_1, title = "Parabacteroides sp. D3 vs D2") -> PARBA_p_1
volcano.plot(table = PARBA_selection_2, title = "Parabacteroides sp. D5 vs D3") -> PARBA_p_2
volcano.plot(table = PARBA_selection_3, title = "Parabacteroides sp. D5 vs D2") -> PARBA_p_3
#BACA, BACB, BACC, CITR, FUS, KER, MONO, PARBA
#D2_D3, D3_D5

pdf(file = "volcano_plot.pdf", width=18, height=16)
grid.arrange(BACA_p_1,BACA_p_2,BACA_p_3,BACB_p_1,BACB_p_2,BACB_p_3,BACC_p_1,BACC_p_2,BACC_p_3,CITR_p_1,CITR_p_2,CITR_p_3,FUS_p_1,FUS_p_2,FUS_p_3,KER_p_1,KER_p_2,KER_p_3,MONO_p_1,MONO_p_2,MONO_p_3,PARBA_p_1,PARBA_p_2,PARBA_p_3, ncol = 3)

dev.off()

pdf(file = "volcano_plot_mono.pdf", width=18, height=16)
grid.arrange(MONO_p_1,MONO_p_2,MONO_p_3, ncol = 3)

dev.off()

pdf(file = "volcano_plot_bacteria.pdf", width=18, height=16)
grid.arrange(BACA_p_1,BACA_p_2,BACA_p_3,BACB_p_1,BACB_p_2,BACB_p_3,BACC_p_1,BACC_p_2,BACC_p_3,CITR_p_1,CITR_p_2,CITR_p_3,FUS_p_1,FUS_p_2,FUS_p_3,KER_p_1,KER_p_2,KER_p_3,PARBA_p_1,PARBA_p_2,PARBA_p_3, ncol = 3)

dev.off()


#EnhancedVolcano(master.selection,
#  lab = master.selection$gene,
#  x='log2FoldChange',
#  y='p.adjusted',
#  title = "Mono D3 vs D2",
#  subtitle = "",
#  pCutoff = 0.05,
#  FCcutoff = 1.0,
#  pointSize = 3.0,
#  legendPosition = 'right',
#  selectLab = c(ifelse(abs(master.selection$log2FoldChange) > 2.5, as.character(master.selection$gene), "")),
#  selectLab = c(ifelse(-log10(master.selection$p.adjusted) > 100 , as.character(master.selection$gene), "")),
#  labSize = 3.0,
#  drawConnectors = TRUE,
#  legendLabSize = 8,
#  legendIconSize = 3,
#  col = c("grey30", "forestgreen", "royalblue", "red2"),
#)
