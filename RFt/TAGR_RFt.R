#TAGR Random Forest model

#https://github.com/Lujun995/RFtest-Simulations.git

#initialization
setwd("RFt/")#to modify according to computer
rm(list=ls())
source("randomForestTest.R")
source("getDescendants.R")
library(GUniFrac)
library(ranger)
library(ape)
library(vegan)
library(parallel)
library(ecodist)
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(randomForest)
library(caret)

set.seed(123)

#loading in data
rooted_tree <- ape::read.tree('tree.nwk')
plot(rooted_tree, node.depth = 5, show.node.label = FALSE, show.tip.label = FALSE)
is.rooted(rooted_tree)

#loading otu table and removing FCSP05 since it has no TTX data
otu_imp <- as.data.frame(read_qza('all_samples_table.qza')$data)
t_otu <- as.data.frame(t(otu_imp)) %>% 
  select(where(is.numeric)) %>% 
  rownames_to_column('SampleID') %>% 
  filter(!SampleID == 'FCSP-05') %>%
  arrange(SampleID) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  .[colSums(.[]) !=0]

#loading metadata
meta.data <- read.csv('newt_metadata.tsv', sep = '\t')
meta.data <- meta.data %>% 
  rename(SampleID = 'X.SampleID') %>% 
  rename(Location = 'Site') %>% 
  filter(!SampleID == 'FCSP-05') %>% 
  arrange(SampleID)

#matching the datasets
filt.meta <- meta.data %>% 
  filter(SampleID %in% rownames(t_otu))

filt.otu <- t_otu %>% 
  rownames_to_column('SampleID') %>% 
  filter(SampleID %in% filt.meta$SampleID) %>%
  arrange(match(SampleID, filt.meta$SampleID))

if (identical(filt.otu$SampleID, filt.meta$SampleID)) {
  print("TRUE")
} else {
  print("FALSE")
}

#join datasets
meta2.data <- filt.meta[, c('SampleID', 'TTX_mg', "Location")]
meta2.data$Location <- as.factor(meta2.data$Location)
meta2.data$Location <- as.integer(meta2.data$Location)

joint.data <- left_join(meta2.data, filt.otu, by = 'SampleID')

joint.data <- joint.data %>% 
  column_to_rownames("SampleID")
joint.data$TTX_mg <- as.numeric(joint.data$TTX_mg)

colnames(joint.data)[-c(1,2)] <- paste0('otu.', colnames(joint.data)[-c(1,2)])

#training the model
tg_x <- joint.data %>% 
  select(!TTX_mg)
y_tg <- joint.data$TTX_mg
index <- createDataPartition(y_tg, p = 0.75, list = FALSE)  

x_train <- tg_x[index, ]
x_test <- tg_x[-index, ]
y_train <- y_tg[index]
y_test <- y_tg[-index]

#actual models
tg_regr <- randomForest(x = x_train, y = y_train, maxnodes = 10, ntree = 100)
predic <- predict(tg_regr, x_test)
result <- x_test
result['TTX_mg'] <- y_test
result['prediction'] <- predic
head(result)

#plotting
plot_data <- data.frame(
  Location = x_test$Location,
  Real = y_test,
  Predicted = predic,
  Loc_actual = meta2.data$Location
)

ggplot() + 
  geom_point(aes(x_test$Location, y = y_test, color = 'red', alpha = 0.5)) +
  geom_point(aes(x_test$Location, y = predic, color = 'blue', alpha = 0.5)) +
  labs(x = 'Location',
       y = 'TTX (mg)',
       color = '') +
  scale_color_manual(labels = c( "Predicted", "Real"), 
                     values = c("blue", "red")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line()
  )
 

plot_data <- data.frame(
  Location = x_test$Location,
  Real = y_test,
  Predicted = predic
)

# Reshape the data for easier plotting
library(tidyr)
plot_data_long <- plot_data %>%
  pivot_longer(cols = c(Real, Predicted), names_to = "Type", values_to = "TTX")

# Create the boxplot with facet wrap 
ggplot(plot_data_long) + 
  geom_boxplot(aes(x = as.factor(Location), y = TTX, color = Type)) +
  labs(x = 'Location',
       y = 'TTX (mg)',
       color = '') +
  scale_color_manual(values = c("Predicted" = "blue", "Real" = "red")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line()
  )

