#initial commit
#LDM
library(usethis)
usethis::use_git_config(user.name = 'Talon-Jost', user.email = 'tallyjost@gmail.com')
rm(list = ls())
#libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(LDM)
library(qiime2R)
library(dplyr)
library(tibble)
# library(ape)
library(mgcv)
library(tidyr)

#import data - this is for Taricha torosa direct with a rarefaction depth of 3800 and an old classifier. proof of concept
metadata <- read.csv('Ttorosa_mapping_Alpha_AntiBd_TTX.csv', header = TRUE) #
metadata <- metadata %>% 
  filter(!Location %in% c('Bolinger', 'Crocker')) %>% 
  filter(!TTX_mg > 2)

OTU_table_3800=as.data.frame(read_qza("rarefied3800_otu-table.qza")$data)

OTU_table_rare=as.data.frame(t(OTU_table_3800))%>%
  rownames_to_column("SampleID")

OTU_table_rare <- as.data.frame(read_qza('rarefied3800_otu-table.qza')$data)

OTU_table<-OTU_table_3800%>%
  rownames_to_column("featureid")

OTU_table_3800_TTX = OTU_table_rare %>%
  filter(SampleID %in% metadata$SampleID)%>%
  column_to_rownames("SampleID")%>%
  .[colSums(.[]) !=0]%>%
  rownames_to_column("SampleID")

OTU_table<-OTU_table_3800%>%
  rownames_to_column("featureid")

OTU_table_filt=OTU_table_3800_TTX%>%
  .[ order(match(OTU_table_3800_TTX$SampleID, metadata$SampleID)), ] %>%
  select(.,-SampleID)

Meta_2 <- metadata %>% 
  filter(SampleID %in% OTU_table_3800_TTX$SampleID)

#taxonomy ####
Taxa<-parse_taxonomy(read_qza("Taxonomy_250.qza")$data)%>%
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

#analysis ####
TTX_ldm <- ldm(OTU_table_filt~Meta_2$Location + Meta_2$TTX_mg,
               perm.within.type="free", 
               perm.between.type="none",
               seed=999,dist.method = "bray", 
               fdr.nominal = 0.01, 
               n.cores = 4)

TTX_ldm$global.tests.stopped
TTX_ldm$otu.tests.stopped
TTX_ldm$p.global.omni
TTX_ldm$detected.otu.omni

DA_OTU_ttx_loc_VE=as.data.frame(TTX_ldm$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
DA_OTU_ttx_loc_F=as.data.frame(TTX_ldm$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
DA_OTU_ttx_loc_qvalue=as.data.frame(TTX_ldm$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_Location=left_join(DA_OTU_ttx_loc_VE,DA_OTU_ttx_loc_F,by="featureid")%>%
  left_join(.,DA_OTU_ttx_loc_qvalue,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxa,by="featureid")

DA_OTU_ttx_TTXp_VE=as.data.frame(DA_OTU_ttx_paper$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
DA_OTU_ttx_TTXp_F=as.data.frame(DA_OTU_ttx_paper$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
DA_OTU_ttx_TTXp_qvalue=as.data.frame(DA_OTU_ttx_paper$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_TTXrank=left_join(DA_OTU_ttx_TTXp_VE,DA_OTU_ttx_TTXp_F,by="featureid")%>%
  left_join(.,DA_OTU_ttx_TTXp_qvalue,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxa,by="featureid") #no sig




#model 2
TTX_ldm2 <- ldm(OTU_table_filt ~ Meta_2$TTX_mg,
                perm.within.type="free", 
                perm.between.type="none",
                seed=999,dist.method = "bray", 
                fdr.nominal = 0.01, 
                n.cores = 4)

TTX_ldm2$global.tests.stopped
TTX_ldm2$otu.tests.stopped
TTX_ldm2$p.global.omni
TTX_ldm2$detected.otu.omni

DA_OTU_ttx_VE2=as.data.frame(TTX_ldm2$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
DA_OTU_ttx_F2=as.data.frame(TTX_ldm2$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
DA_OTU_ttx_qvalue2=as.data.frame(TTX_ldm2$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_TTX=left_join(DA_OTU_ttx_VE2,DA_OTU_ttx_F2,by="featureid")%>%
  left_join(.,DA_OTU_ttx_qvalue2,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxa,by="featureid") #no sig ###

#model 3 - just location
TTX_ldm3 <- ldm(OTU_table_filt ~ Meta_2$Location,
                perm.within.type="free", 
                perm.between.type="none",
                seed=999,dist.method = "bray", 
                fdr.nominal = 0.01, 
                n.cores = 4)

TTX_ldm3$global.tests.stopped
TTX_ldm3$otu.tests.stopped
TTX_ldm3$p.global.omni
TTX_ldm3$detected.otu.omni

DA_OTU_loc_VE3=as.data.frame(TTX_ldm3$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
DA_OTU_loc_F3=as.data.frame(TTX_ldm3$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
DA_OTU_loc_qvalue3=as.data.frame(TTX_ldm3$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_loc3=left_join(DA_OTU_loc_VE3,DA_OTU_loc_F3,by="featureid")%>%
  left_join(.,DA_OTU_loc_qvalue3,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxa,by="featureid") #there apears to be about 30 significant OTUs


#model 4 - location as significant taking into account TTX concentration and Bd infection status
TTX_ldm4 <- ldm(OTU_table_filt | (Meta_2$TTX_mg + Meta_2$Infected) ~ Meta_2$Location,
                perm.within.type="free", 
                perm.between.type="none",
                seed=999,dist.method = "bray", 
                fdr.nominal = 0.01, 
                n.cores = 4)

TTX_ldm4$global.tests.stopped
TTX_ldm4$otu.tests.stopped
TTX_ldm4$p.global.omni
TTX_ldm4$detected.otu.omni

DA_OTU_loc_ttx_sex_VE4=as.data.frame(TTX_ldm4$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
DA_OTU_loc_ttx_sex_F4=as.data.frame(TTX_ldm4$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
DA_OTU_loc_ttx_sex_qvalue4=as.data.frame(TTX_ldm4$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_loc3=left_join(DA_OTU_loc_ttx_sex_VE4,DA_OTU_loc_ttx_sex_F4,by="featureid")%>%
  left_join(.,DA_OTU_loc_ttx_sex_qvalue4,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxa,by="featureid") #there apears to be about 15 significant OTUs when accounting for Bd infection and TTX concentration for location
