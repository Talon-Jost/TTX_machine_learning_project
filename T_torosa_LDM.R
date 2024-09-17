#initial commit
#LDM

rm(list = ls())
#libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(LDM)
library(qiime2R)
library(dplyr)
library(tibble)
library(ape)
library(mgcv)
library(tidyr)

#import data
metadata <- read.csv("tt_metadata.tsv", sep = '\t') %>% 
  as.data.frame()%>%
  filter(TTX_Richness != is.na(TTX_Richness))

metadata_TTX=metadata%>%
  filter(!Location%in%c("Bolinger","Crocker"),
         !TTX_mg>2)

OTU_table_2200=as.data.frame(read_qza("rarefied_table.qza")$data)

OTU_table_rare=as.data.frame(t(OTU_table_2200))%>%
  rownames_to_column("SampleID")

#Taxonomy
Taxa<-parse_taxonomy(read_qza("classification.qza")$data) %>% 
  rownames_to_column("featureid") %>% 
  dplyr::select(-Kingdom)%>%
  mutate(across(everything(),~gsub("[][]","",.)),
         across(everything(),~replace_na(.,"")))

for (i in 1:nrow(Taxa)){
  if (Taxa[i,7] != ""){
    Taxa$Species[i] <- paste(Taxa$Genus[i], Taxa$Species[i], sep=" ")
  } else if (Taxa[i,2] == ""){
    Phylum <- paste("Unclassified", "Bacteria", sep="-")
    Taxa[i, 2:7] <- Phylum
  } else if (Taxa[i,3] == ""){
    Class <- paste(Taxa[i,2],"(P)", sep=" ")
    Taxa[i, 3:7] <- Class
  } else if (Taxa[i,4] == ""){
    Order <- paste(Taxa[i,3],"(C)", sep=" ")
    Taxa[i, 4:7] <- Order
  } else if (Taxa[i,5] == ""){
    Family <- paste(Taxa[i,4],"(O)", sep=" ")
    Taxa[i, 5:7] <- Family
  } else if (Taxa[i,6] == ""){
    Genus <- paste(Taxa[i,5],"(F)", sep=" ")
    Taxa[i, 6:7] <- Genus
  } else if (Taxa[i,7] == ""){
    Taxa$Species[i] <- paste(Taxa$Genus[i],"sp.", sep=" ")
  }
}

# Differential Abundance Analysis using LDM ####
#filter OTU table to metadata
OTU_table_2200_TTX = OTU_table_rare %>%
  filter(SampleID %in% metadata_TTX$SampleID)%>%
  column_to_rownames("SampleID")%>%
  .[colSums(.[]) !=0]%>%
  rownames_to_column("SampleID")

#match, filter, and order metadata and OTU table
OTU_table_filt=OTU_table_2200_TTX %>%
  .[ order(match(OTU_table_2200_TTX$SampleID, metadata_TTX$SampleID)), ] %>% 
  as.data.frame()
  
rownames(OTU_table_filt) = OTU_table_filt$SampleID
OTU_table_filt = OTU_table_filt[,-1]

# check_tab1=OTU_table_2200%>%
#   rownames_to_column(var = 'featureid')
# check_tab2 <- check_tab1 %>% 
#   .[ order(match(check_tab1$featureid, Taxa$featureid)), ]

Meta_2 <- metadata_TTX %>% 
  filter(SampleID %in% OTU_table_rare$SampleID)

#set up model
TTX_ldm <- ldm(OTU_table_filt~Meta_2$Location + Meta_2$TTX_mg,
               perm.within.type="free", perm.between.type="none",
               seed=999,dist.method = "bray", fdr.nominal = 0.01, n.cores = 4)

TTX_ldm$global.tests.stopped
TTX_ldm$otu.tests.stopped
TTX_ldm$p.global.omni
TTX_ldm$detected.otu.omni


r1 <- match(TTX_ldm$detected.otu.omni[[1]], colnames(TTX_ldm$q.otu.omni))
order1 <- order(TTX_ldm$p.otu.omni[1,])[1:10]
order1 <- as.data.frame(order1)

sum.tab <- data.frame(
  raw.pvalue = signif(TTX_ldm$p.otu.omni[1, order1], 5),
  adj.pvalue = signif(TTX_ldm$q.otu.omni[1, order1], 5),
  mean.freq = signif(TTX_ldm$mean.freq[order1], 5),
  direction = t(ifelse(TTX_ldm$beta[1,] > 0, "+", "-"))[order1,],
  featureid = colnames(TTX_ldm$q.otu.omni)[order1],
  row.names = NULL)

print(sum.tab)
#no significance with this model with the adjusted p-val so we will leave this for later

LDM_TTXrank1=left_join(TTX_ldm_ve,ttx_ldm_p_f,by="featureid")%>%
  left_join(.,TTX_ldm_q_value,by="featureid")%>%
  filter(!q_value_cov2 > 0.050)%>%
  left_join(.,Taxa,by="featureid")






#Model 2
TTX_ldm2 <- ldm(OTU_table_filt~ Meta_2$TTX_mg,
               perm.within.type="free", perm.between.type="none",
               seed=999,dist.method = "bray", fdr.nominal = 0.01, n.cores = 4)

TTX_ldm2$global.tests.stopped
TTX_ldm2$otu.tests.stopped
TTX_ldm2$p.global.omni
TTX_ldm2$detected.otu.omni
TTX_ldm2$p.otu.omni
TTX_ldm2$q.otu.omni
dim(TTX_ldm2$beta)

TTX_ldm2$q.otu.omni[is.na(TTX_ldm2$q.otu.omni)] <- 0
TTX_ldm2$mean.freq[is.na(TTX_ldm2$mean.freq)] <- 0
TTX_ldm2$p.otu.omni[is.na(TTX_ldm2$p.otu.omni)] <- 0

r2 <- match(TTX_ldm2$detected.otu.omni[[1]], colnames(TTX_ldm2$q.otu.omni))
order2 <- as.numeric(order(TTX_ldm2$p.otu.omni)[1])
order2 <- as.numeric(as.data.frame(order2))

pvalues <- TTX_ldm2$p.otu.omni[[1]]

# Ensure order is numeric
order2 <- as.numeric(order(pvalues)[1:10])


sum.tab2 <- data.frame(
  raw.pvalue = signif(pvalues[order2], 10),
  adj.pvalue = signif(TTX_ldm2$q.otu.omni[,order2], 10),
  mean.freq = signif(TTX_ldm2$mean.freq[,order2], 10),
  direction = t(ifelse(TTX_ldm2$beta> 0, "+", "-"))[order2],
  featureid = names(TTX_ldm2$q.otu.omni)[order2],
  row.names = NULL)

print(sum.tab2)





#Test 3 - just location

TTX_ldm3 <- ldm(OTU_table_filt~ Meta_2$Location,
                perm.within.type="free", perm.between.type="none",
                seed=999,dist.method = "bray", fdr.nominal = 0.01, n.cores = 4)

TTX_ldm3$global.tests.stopped
TTX_ldm3$otu.tests.stopped
TTX_ldm3$p.global.omni
TTX_ldm3$detected.otu.omni
TTX_ldm3$p.otu.omni
TTX_ldm3$q.otu.omni
dim(TTX_ldm3$beta)

TTX_ldm3$q.otu.omni[is.na(TTX_ldm3$q.otu.omni)] <- 0
TTX_ldm3$mean.freq[is.na(TTX_ldm3$mean.freq)] <- 0
TTX_ldm3$p.otu.omni[is.na(TTX_ldm3$p.otu.omni)] <- 0

r3 <- match(TTX_ldm3$detected.otu.omni[[1]], colnames(TTX_ldm3$q.otu.omni))
order3 <- as.numeric(order(TTX_ldm3$p.otu.omni)[1])
order3 <- as.numeric(as.data.frame(order3))

pvalues2 <- TTX_ldm3$p.otu.omni[[1]]

# Ensure order is numeric
order3 <- as.numeric(order(pvalues2)[1:10])


sum.tab3 <- data.frame(
  raw.pvalue = signif(pvalues2[order2], 10),
  adj.pvalue = signif(TTX_ldm3$q.otu.omni[,order3], 10),
  mean.freq = signif(TTX_ldm3$mean.freq[,order3], 10),
  direction = t(ifelse(TTX_ldm3$beta> 0, "+", "-"))[order3],
  featureid = names(TTX_ldm3$q.otu.omni)[order3],
  row.names = NULL)

print(sum.tab3)




#check to see if they're different

BC_LDM1 <- permanovaFL(tg_bray_curtis_matrix ~ filtered_metadata$Location, 
                       perm.between.type = 'free',
                       perm.within.type = 'free',
                       seed = 999,
                       square.dist = FALSE,
                       n.cores = 8,
                       n.perm.max = 5000)
BC_LDM1$F.statistics
BC_LDM1$R.squared
BC_LDM1$p.permanova

#how are they different

BC_disp<-betadisper(as.dist(tg_bray_curtis_matrix), 
                    filtered_metadata$Location,
                    type = c("centroid"),
                    bias.adjust = TRUE,
                    sqrt.dist = FALSE,
                    add = FALSE)
anova(BC_disp)


Taxa <-parse_taxonomy(read_qza("classification.qza")$data)%>%
  rownames_to_column("featureid")%>%
  dplyr::select(-Kingdom)%>%
  mutate(across(everything(),~gsub("[][]","",.)),
         across(everything(),~replace_na(.,"")))

for (i in 1:nrow(Taxa)){
  if (Taxa[i,7] != ""){
    Taxa$Species[i] <- paste(Taxa$Genus[i], Taxa$Species[i], sep=" ")
  } else if (Taxa[i,2] == ""){
    Phylum <- paste("Unclassified", "Bacteria", sep="-")
    Taxa[i, 2:7] <- Phylum
  } else if (Taxa[i,3] == ""){
    Class <- paste(Taxa[i,2],"(P)", sep=" ")
    Taxa[i, 3:7] <- Class
  } else if (Taxa[i,4] == ""){
    Order <- paste(Taxa[i,3],"(C)", sep=" ")
    Taxa[i, 4:7] <- Order
  } else if (Taxa[i,5] == ""){
    Family <- paste(Taxa[i,4],"(O)", sep=" ")
    Taxa[i, 5:7] <- Family
  } else if (Taxa[i,6] == ""){
    Genus <- paste(Taxa[i,5],"(F)", sep=" ")
    Taxa[i, 6:7] <- Genus
  } else if (Taxa[i,7] == ""){
    Taxa$Species[i] <- paste(Taxa$Genus[i],"sp.", sep=" ")
  }
}

rar_tab <- as.data.frame(read_qza('./core_metrics/all_samples/2200/all_samples_table.qza')$data)
tg_2200_rare <- as.data.frame(t(rar_tab)) %>%
  rownames_to_column("SampleID")

tg_OTU_table <- rar_tab %>%
  rownames_to_column('featureID')

tg_2200_rare <- tg_2200_rare %>%
  filter(!SampleID == 'FCSP-05') %>%
  arrange("SampleID")

matched_1 <- tg_2200_rare %>%
  filter(SampleID %in% filtered_metadata$SampleID)%>%
  column_to_rownames("SampleID")%>%
  .[colSums(.[]) !=0] %>%
  rownames_to_column("SampleID")

matched_tab <- matched_1 %>%
  .[ order(match(matched_1$SampleID, metadata_frame$SampleID)),] %>%
  select(.,-SampleID)

metadata_frame <- read.csv('~/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/newt_metadata.tsv', 
                           sep = '\t',
                           row.names = NULL) %>%
  rename(SampleID = 1) %>%
  filter(!SampleID == 'FCSP-05') %>%
  as.data.frame() %>%
  arrange(SampleID) %>%
  filter(SampleID %in% filtered_metadata$SampleID)

metadata_new <- metadata_frame %>%
  select(SampleID, TTX_mg) %>%
  column_to_rownames("SampleID")

ttx_mg <- as.data.frame(metadata_new$TTX_mg)
colnames(ttx_mg)[1] <- "TTX_mg"
rownames_to_column(ttx_mg)
rownames(ttx_mg) <- NULL
any(duplicated(rownames(ttx_mg)))
any(duplicated(rownames(matched_tab)))
rownames(metadata_frame) <- seq_len(nrow(metadata_frame))

tg_ttx_location_ldm <- ldm(matched_tab ~ metadata_frame$Site,
                           perm.within.type="free", 
                           perm.between.type="none",
                           seed=999,dist.method = "bray", 
                           fdr.nominal = 0.01, 
                           n.cores = 10,
                           n.perm.max = 10000)