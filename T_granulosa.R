
#Taricha granulosa analysis using LDM
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

#Metadata importation. FCSP-05 doesn't have a TTX measurement so it will be removed entirely for this analysis
tagr_metadata <- read.csv('newt_metadata.tsv', sep = '\t', row.names = NULL)
meta <- tagr_metadata %>%
  filter(!X.SampleID %in% c('Pos', 'Neg', 'FCSP-05')) %>% 
  rename(SampleID = 'X.SampleID') %>%
  rename(Location = "Site") %>% 
  arrange(SampleID)

meta$TTX_mg <- as.numeric(meta$TTX_mg)
meta$Location <- as_factor(meta$TTX_mg)

#OTU table
rarefied_otu_table <- as.data.frame(read_qza('all_samples_table.qza')$data)
otu_table <- as.data.frame(t(rarefied_otu_table)) %>% 
  select(where(is.numeric)) %>%
  rownames_to_column("SampleID") %>% 
  filter(!SampleID %in% "FCSP-05") %>% 
  column_to_rownames(var = "SampleID") %>% 
  .[colSums(.[]) !=0]%>% 
  rownames_to_column("SampleID")

#filtering those to match with confirmation
filt_metadata <- meta %>% 
  filter(SampleID %in% otu_table$SampleID)

filt_otu_table2 <- otu_table %>% 
  filter(SampleID %in% filt_metadata$SampleID) %>%
  arrange(match(SampleID, filt_metadata$SampleID))

if (identical(filt_otu_table2$SampleID, filt_metadata$SampleID)) {
  print("TRUE")
} else {
  print("FALSE")
}

rownames(filt_otu_table2) <- NULL
rownames(filt_metadata) <- NULL
filt_otu_table2 <- filt_otu_table2 %>% 
  select(-SampleID)
filt_metadata[40, "TTX_mg"] <- filt_metadata[40, "TTX_mg"] + 0.0001
filt_metadata[64, "TTX_mg"] <- filt_metadata[64, "TTX_mg"] + 0.0001
filt_metadata[71, "TTX_mg"] <- filt_metadata[71, "TTX_mg"] + 0.0001
filt_metadata[66, "TTX_mg"] <- filt_metadata[66, "TTX_mg"] + 0.0001
filt_metadata[81, "TTX_mg"] <- filt_metadata[81, "TTX_mg"] + 0.0001
#Taxonomy ###
Taxa <- as.data.frame(read_qza("classification.qza")$data) %>% 
  select(!Confidence) %>% 
  mutate(across(everything(),~gsub("[][]","",.)),
         across(everything(),~replace_na(.,"")))

Taxonomy <- Taxa %>%
  mutate(taxonomy_parts = str_split(Taxon, "; ")) %>%
  unnest(taxonomy_parts) %>%
  separate(taxonomy_parts, into = c('Rank', 'Value'), sep = "__") %>%
  pivot_wider(names_from = Rank, values_from = Value) %>% 
  rename(
    Domain = d, Phylum = p, Class = c,
    Order = o, Family = f, Genus = g,
    Species = s
  ) %>% 
  select(-Taxon) %>% 
  column_to_rownames(var = "Feature.ID")

Taxonomy <- Taxonomy %>% 
  mutate(Species = if_else(is.na(Species), 'sp.', Species)) %>% 
  mutate(Genus = if_else(is.na(Genus), paste(Family, 'Genus'), Genus)) %>% 
  rownames_to_column("featureid")

#LDM1 model with TTX and location
tagr_TTX_ldm <- ldm(filt_otu_table2 ~ filt_metadata$TTX_mg,
               perm.within.type="free", 
               perm.between.type="none",
               seed=999,dist.method = "bray", 
               fdr.nominal = 0.01, 
               n.cores = 4)

tagr_TTX_ldm$global.tests.stopped
tagr_TTX_ldm$otu.tests.stopped
tagr_TTX_ldm$p.global.omni
tagr_TTX_ldm$detected.otu.omni

tagr_DA_OTU_ttx_loc_VE=as.data.frame(tagr_TTX_ldm$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
tagr_DA_OTU_ttx_loc_F=as.data.frame(tagr_TTX_ldm$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
tagr_DA_OTU_ttx_loc_qvalue=as.data.frame(tagr_TTX_ldm$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_Location=left_join(tagr_DA_OTU_ttx_loc_VE,tagr_DA_OTU_ttx_loc_F,by="featureid")%>%
  left_join(.,tagr_DA_OTU_ttx_loc_qvalue,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxonomy,by="featureid")


LDM_TTXrank=left_join(DA_OTU_ttx_TTXp_VE,DA_OTU_ttx_TTXp_F,by="featureid")%>%
  left_join(.,DA_OTU_ttx_TTXp_qvalue,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxonomy,by="featureid") #11 sig values

  
#LDM location ####
tagr_TTX_ldm <- ldm(filt_otu_table2 ~ filt_metadata$Location,
                    perm.within.type="free", 
                    perm.between.type="none",
                    seed=999,dist.method = "bray", 
                    fdr.nominal = 0.01, 
                    n.cores = 4)

tagr_TTX_ldm$global.tests.stopped
tagr_TTX_ldm$otu.tests.stopped
tagr_TTX_ldm$p.global.omni
tagr_TTX_ldm$detected.otu.omni

tagr_DA_OTU_ttx_loc_VE=as.data.frame(tagr_TTX_ldm$VE.otu.freq.submodels)%>%
  rownames_to_column("featureid")%>%
  rename(VE=names(.)[2])
tagr_DA_OTU_ttx_loc_F=as.data.frame(tagr_TTX_ldm$F.otu.freq)%>%
  rownames_to_column("featureid")%>%
  rename(F_value=names(.)[2])
tagr_DA_OTU_ttx_loc_qvalue=as.data.frame(tagr_TTX_ldm$q.otu.omni)%>%
  rownames_to_column("featureid")%>%
  rename(q_value=names(.)[2])

LDM_Location=left_join(tagr_DA_OTU_ttx_loc_VE,tagr_DA_OTU_ttx_loc_F,by="featureid")%>%
  left_join(.,tagr_DA_OTU_ttx_loc_qvalue,by="featureid")%>%
  filter(!q_value > 0.050)%>%
  left_join(.,Taxonomy,by="featureid")