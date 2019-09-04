#!/usr/bin/env Rscript
#argv <- commandArgs(TRUE)
#z <- c(argv[1])
#setwd(dir = z)
setwd("C:/Users/admin/Desktop/TF/data")
essential.data <- read.csv("db_nrg3225.txt",sep ="\t",header = T)
essential_GID <- essential.data$G_ID
HIPred.data <- read.csv("HIPred_IDset.csv",sep ="\t",header = F)
HIPred_GID <- HIPred.data$V3
nature.data <- read.csv("db_2nature3230.txt",sep ="\t",header = T)
nature_GID <- nature.data$transcript
CCN.data <- read.csv("db_data7014.txt",sep ="\t",header = T)
CCN_GID <- CCN.data$Ensembl.Gene.ID
Ohnolog.data <- read.csv("db_data7294.txt",sep ="\t",header = T)
Ohnolog_GID <- Ohnolog.data$Ensembl.id
title <- c("essential_GID",
           "HIPred_GID",
           "nature_GID",
           "CCN_GID",
           "Ohnolog_GID"
)
#nuclear_genome
gene.data <- read.csv("dataset_nonTID.txt",sep =",",header = T)
chr <- c(1:22,"X","Y")
human_chr_GID <- gene.data[gene.data$Chromosome.scaffold.name %in% chr, ]$Gene.stable.ID
human.data <- gene.data[gene.data$Chromosome.scaffold.name %in% chr, ]
#TF_in_gene
trans.data <- read.csv("TableS1.csv",sep =",",header = T)
TF_GID <- trans.data[trans.data$Is.TF. == "Yes",]$Gene.Information
TF.data <- trans.data[trans.data$Is.TF. == "Yes",]
#split_human_TF_nTF
nTF <- setdiff(human_chr_GID, TF_GID)
TF <- intersect(human_chr_GID, TF_GID)
#DBD
DBD <- TF.data$X.1
DBD_family <- unique(DBD)
#hash_DBD_TF_GID
#install.packages("hash")
library(hash)
GID_TF_DBD <- hash()
#DBD_TF_GID <- hash()
values <- NULL
for (Gene.Information in TF.data) {
  
  .set(GID_TF_DBD, keys = TF.data$Gene.Information, values = TF.data$X.1)
}
#DBD_TF_GID <- rev(GID_TF_DBD)
#
#fiher_test_factor_DBD
#DBD.data <- with(TF.data, data.frame(GID = Gene.Information,
#                                     Gname = x,
#                                     DBD = x.1))


#Fisher_test_factor_all
for (sample in title){
  sample_TF <- intersect(TF, sample)
  non_sample_TF <- setdiff(TF, sample)
  sample_nTF <- intersect(nTF, sample)
  non_sample_nTF <- setdiff(nTF, sample)
  fisher_sample <- matrix(data = c(length(sample_TF),
                                       length(non_sample_TF),
                                       length(sample_nTF),
                                       length(non_sample_nTF), nrow = 2))
  colnames(fisher_sample) <- c("TF", "non_TF")
  rownames(fisher_sample) <- c("Factor", "non_Factor")
  fisher.test(fisher_sample)
}


#example£ºFisher_test_factor_all
Ohnolog_TF <- intersect(TF, Ohnolog_GID)
non_Ohnolog_TF <- setdiff(TF, Ohnolog_GID)
Ohnolog_nTF <- intersect(nTF, Ohnolog_GID)
non_Ohnolog_nTF <- setdiff(nTF, Ohnolog_GID)
fisher_Ohnolog <-matrix(data = c(length(Ohnolog_TF),
                             length(non_Ohnolog_TF),
                             length(Ohnolog_nTF),
                             length(non_Ohnolog_nTF)),nrow = 2)
colnames(fisher_Ohnolog) <- c("TF", "non_TF")
rownames(fisher_Ohnolog) <- c("Factor", "non_Factor")
fisher.test(fisher_Ohnolog)
result <- fisher.test(fisher_Ohnolog)



#TF_five_times
TF_essensial <- intersect(TF_GID, essential_GID)
TF_HIPred <- intersect(TF_GID, HIPred_GID)
TF_nature <- intersect(TF_GID, nature_GID)
TF_CCN <- intersect(TF_GID, CCN_GID)
TF_Ohnolog <- intersect(TF_GID, Ohnolog_GID)

TF_twice <- intersect(TF_essensial, TF_HIPred)
TF_third <- intersect(TF_twice, TF_nature)
TF_fourth <- intersect(TF_third, TF_CCN)
TF_fifth <- intersect(TF_fourth, TF_Ohnolog)
write.csv(TF_fifth, "TF_fifth.csv", row.names = FALSE, quote = F)
#nTF_five_times
TF_non_Ohnolog <- setdiff(TF_GID, Ohnolog_GID)
TF_non_essential <- setdiff(TF_GID, essential_GID)
TF_non_HIPred <- setdiff(TF_GID, HIPred_GID)
TF_non_nature <- setdiff(TF_GID, nature_GID)
TF_non_CCN <- setdiff(TF_GID, CCN_GID)
TF_non_twice <- intersect(TF_non_Ohnolog, TF_non_essential)
TF_non_third <- intersect(TF_non_twice, TF_non_HIPred)
TF_non_fourth <- intersect(TF_non_third, TF_non_nature)
TF_non_fifth <- intersect(TF_non_fourth, TF_non_CCN)
write.csv(TF_non_fifth, "TF_non_fifth.csv", row.names = FALSE, quote = F)
#neutral_TF = TF_GID - TF_fifth - TF_non_fifth
nfifth_TF <- setdiff(TF_GID, TF_fifth)
neutral_TF <- setdiff(nfifth_TF, TF_non_fifth)
#dN_dS_with_mouse
with_mouse.data <- read.csv("dN_dS_with_mouse.csv", sep =",", header = T)
#one_copy_number
with_mouse.data <- unique(with_mouse.data)
index <- duplicated(with_mouse.data$Gene.stable.ID)
#dN_dS.data2 <- dN_dS.data[!index,]
index2 <- duplicated(with_mouse.data$Mouse.gene.stable.ID)
#dN_dS.data3 <- dN_dS.data[!index2,]
OCN.data <- with_mouse.data[(!index&!index2),]
write.csv(OCN.data, "one_copy_number.csv", row.names = FALSE, quote = F)
#dN/dS
dN_dS.data <- with(OCN.data, 
      
                                data.frame(Gene.stable.ID = Gene.stable.ID, 
                              Mouse.gene.stable.ID= Mouse.gene.stable.ID,
                              dN_dS = dN.with.Mouse / dS.with.Mouse))
#dN_dS_TF_fifth.data
TF_fifth.dN_dS <- NULL
for (gid in TF_fifth) {
  TF_fifth.dN_dS <- rbind(TF_fifth.dN_dS, dN_dS.data[dN_dS.data$Gene.stable.ID == gid, ])
}
write.csv(TF_fifth.dN_dS, "TF_fifth_dN_dS.csv", row.names = FALSE, quote = F)
#dN_dS_TF_non_fifth.data
TF_non_fifth.dN_dS <- NULL
for (gid in TF_non_fifth) {
  TF_non_fifth.dN_dS <- rbind(TF_non_fifth.dN_dS, dN_dS.data[dN_dS.data$Gene.stable.ID == gid, ])
}
write.csv(TF_non_fifth.dN_dS, "TF_fifth_non_dN_dS.csv", row.names = FALSE, quote = F)
#dN_dS_TF_neutral.data
neutral_TF.dN_dS <- NULL
for (gid in neutral_TF) {
  neutral_TF.dN_dS <- rbind(neutral_TF.dN_dS, dN_dS.data[dN_dS.data$Gene.stable.ID == gid, ])
}
write.csv(neutral_TF.dN_dS, "neutral_TF_dN_dS.csv", row.names = FALSE, quote = F)


