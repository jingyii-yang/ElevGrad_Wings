setwd("X:/home/_Wings_/HPC/")

library(tidyverse)
library(ape)
library(caper)
library(phylolm)


# 100 randomly selected tree from 10,000 full trees (9,993 species; Hackett backbone)
n_trees <- readRDS('100 random trees.rds') 
# read the BirdTree version of dataset
dat = read.csv('Data S1_BirdTree.csv') %>% rename(MAX = 'Max.Elevation.1', MEAN = 'Mean.Elevation.1')
# remove flightless (and NA) species - irrelevant to hypothesis
dat = filter(dat, Flight.Mode != 'flightless')
# align species names to phylogeny
dat$Species_BirdTree <- gsub(' ', '_', dat$Species_BirdTree) 
    


##### Data transformation #####

## normalise and standardise data
{
  dat$log.Mass <- scale(log(dat$Body.Mass))/2
  
  dat$log.Wing.Length <- scale(log(dat$Wing.Length))/2
  dat$log.Sec1 <- scale(log(dat$Secondary1))/2
  dat$log.WA = scale(log(dat$Wing.Length * dat$Secondary1 * pi/2))/2
  dat$HWI <- scale(dat$Hand.Wing.Index)/2
  
  # Environmental continuous variables will be transformed later, within each elevation band, 
  # to make their effect sizes comparable across bands (especially for the main predictor Elevation).
}


## assign elevation bands to each species
{
  dat$elevbin1 = ifelse(dat$MAX <= 3000, 1, 0) 
  dat$elevbin2 = ifelse(dat$MAX >= 1000 & dat$MAX <= 4000, 1, 0)
  dat$elevbin3 = ifelse(dat$MAX >= 2000 & dat$MAX <= 5000, 1, 0)
  dat$elevbin4 = ifelse(dat$MAX >= 3000 & dat$MAX <= 6000, 1, 0)
  dat$elevbin5 = ifelse(dat$MAX >= 4000 & dat$MAX <= 7000, 1, 0)
  dat$elevbin6 = ifelse(dat$MAX >= 5000, 1, 0) 
}


########################################### MODELLING #####################################################


## PART 1: do models within each elevation band, use the same tree for each of the 4 metrics

elev_bin <- function(treei, dataset, elev_type){
  
  # choose elevation band columns
  bin_cols <<- str_which(colnames(dataset), 'elevbin')
  
  coef = c()
  for (i in bin_cols){
    
    ## select species in the current band (i)
    dat_bin = dataset %>% filter_at(i, all_vars(. == 1))
    
    {
      ## re-scale Elev to keep the data spread consistent (i.e. bandwidth is consistent throughout, so should be the predictor Elev), 
      ## then the results will be comparable across bands.
      dat_bin$Elev = scale(sqrt(dat_bin[, elev_type]))/2
      dat_bin$Lat <- scale(sqrt(abs(dat_bin$Centroid.latitude)))/2
      dat_bin$TempVar <- scale(log(dat_bin$Temp..seasonality))/2
      
      ## NO need to re-scale other continuous variables because:
      # (A) the data spread of wing metrics (Ys), wider or narrower in the current band, reflects reality
      #     (e.g. if HWI has lower spread in 3-6km it means HWI DOES have less variation there than in other elevation bands);
      # (B) globally all Ys have been scaled so the results are still comparable across wing metrics;
      # (C) do not re-scale body mass because each mass data point pairs with its wing area/length/width, so it must be treated the same way as Ys.
    }
    
  
    datree_bin <- caper::comparative.data(data = dat_bin, phy = treei,          
                                          names.col = 'Species_BirdTree',
                                          na.omit = FALSE)
    
    
    res_hwi <- phylolm(HWI ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                       data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_hwi)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_hwi$formula[[2]])
    ## sample size depends on flight mode because adding the interaction means 
    ## the slope estimate can represent for the reference group only. 
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]  
    coef = rbind(coef, coef_i)
    
    
    res_wa <- phylolm(log.WA ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wa)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wa$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
    
    
    res_wl <- phylolm(log.Wing.Length ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
    
    
    res_sl <- phylolm(log.Sec1 ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_sl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_sl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
    
  }
  return(coef)
}



## PART 2: loop through the 100 trees and merge all results

elev_bin_100coef <- function(input, tag, input_elev = 'MAX'){
  
  coef_ntree = c()
  for (ii in 1:100){
    coef_treeii = elev_bin(treei = n_trees[[ii]], dataset = input, elev_type = input_elev)
    coef_treeii$treeID = ii
    # this table will contain all coefficients from all elevation bands x 4 wing metrics x 100 trees
    coef_ntree = rbind(coef_ntree, coef_treeii)
  }
    saveRDS(coef_ntree, paste0('X:/home/_Wings_/HPC/elev_band_coef100/', tag, '.rds'))
  return(coef_ntree)
}



########  Main analysis [Figure 4] #########

elev_bin_100coef(input = dat, tag = 'all species')


######  Supplementary analysis: use bandwidth = 2km [Figure S4B] ##############

dat_2km = dat

## re-assign elevation bands to each species
{ # Remove existing 'elevbin' columns
  dat_2km = dplyr::select(dat_2km, -contains('elevbin'))
  
  dat_2km$elevbin1 = ifelse(dat_2km$MAX <= 2000, 1, 0) 
  dat_2km$elevbin2 = ifelse(dat_2km$MAX >= 1000 & dat_2km$MAX <= 3000, 1, 0)
  dat_2km$elevbin3 = ifelse(dat_2km$MAX >= 2000 & dat_2km$MAX <= 4000, 1, 0)
  dat_2km$elevbin4 = ifelse(dat_2km$MAX >= 3000 & dat_2km$MAX <= 5000, 1, 0)
  dat_2km$elevbin5 = ifelse(dat_2km$MAX >= 4000 & dat_2km$MAX <= 6000, 1, 0)
  dat_2km$elevbin6 = ifelse(dat_2km$MAX >= 5000 & dat_2km$MAX <= 7000, 1, 0)
}
elev_bin_100coef(input = dat_2km, tag = 'bandwidth = 2km')




######  Supplementary analysis: use bandwidth = 4km [Figure S4C]  ###############

dat_4km = dat

## re-assign elevation bands to each species
{ # Remove existing 'elevbin' columns
  dat_4km = dplyr::select(dat_4km, -contains('elevbin'))
  
  dat_4km$elevbin1 = ifelse(dat_4km$MAX <= 4000, 1, 0)
  dat_4km$elevbin2 = ifelse(dat_4km$MAX >= 1000 & dat_4km$MAX <= 5000, 1, 0)
  dat_4km$elevbin3 = ifelse(dat_4km$MAX >= 2000 & dat_4km$MAX <= 6000, 1, 0)
  dat_4km$elevbin4 = ifelse(dat_4km$MAX >= 3000 & dat_4km$MAX <= 7000, 1, 0)
  dat_4km$elevbin5 = ifelse(dat_4km$MAX >= 4000, 1, 0)
}
elev_bin_100coef(input = dat_4km, tag = 'bandwidth = 4km')



######  Supplementary analysis: non-migratory landbirds [[Figure S4D] ] ###############

dat2 <- dat %>% filter(Migration != 3) %>% filter(Habitat.type != 'Marine')

elev_bin_100coef(input = dat2, tag = 'landbird model')



######  Supplementary analysis: retest with high certainty data [Figure S4E]  ######

dat_A = filter(dat, Aerial.Lifestyle.Certainty == 'A') 

elev_bin_100coef(input = dat_A, tag = 'high certainty')



######  Supplementary analysis: retest with mean elevation data [Figure S4F]  ######

dat_mean = dat

## re-assign elevation bands based on mean elevation
{ # Remove existing 'elevbin' columns
  dat_mean = dplyr::select(dat_mean, -contains('elevbin'))
  
  dat_mean$elevbin1 = ifelse(dat_mean$MEAN <= 3000, 1, 0)
  dat_mean$elevbin2 = ifelse(dat_mean$MEAN >= 1000 & dat_mean$MEAN <= 4000, 1, 0)
  dat_mean$elevbin3 = ifelse(dat_mean$MEAN >= 2000 & dat_mean$MEAN <= 5000, 1, 0)
  dat_mean$elevbin4 = ifelse(dat_mean$MEAN >= 3000, 1, 0) 
}

elev_bin_100coef(input = dat_mean, tag = 'mean.elevation', input_elev = 'MEAN')



######  Supplementary analysis: absolute wing metric changes [Figure S4A]  ######

## (models/functions are nearly the same except that body mass needs to be removed from predictors)

## PART 1: do models within each elevation band, use the same tree for each of the 4 metrics
elev_bin_abs <- function(treei, dataset){
  
  bin_cols <<- str_which(colnames(dataset), 'elevbin')
  
  coef = c()
  for (i in bin_cols){
    
    dat_bin = dataset %>% filter_at(i, all_vars(. == 1))
    
    {dat_bin$Elev = scale(sqrt(dat_bin$MAX))/2
     dat_bin$Lat <- scale(sqrt(abs(dat_bin$Centroid.latitude)))/2
     dat_bin$TempVar <- scale(log(dat_bin$Temp..seasonality))/2}
    
    datree_bin <- caper::comparative.data(data = dat_bin, phy = treei,          
                                          names.col = 'Species_BirdTree',
                                          na.omit = FALSE)
    
    res_hwi <- phylolm(HWI ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level,
                       data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_hwi)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_hwi$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]  
    coef = rbind(coef, coef_i)
    
    res_wa <- phylolm(log.WA ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wa)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wa$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
    
    res_wl <- phylolm(log.Wing.Length ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
    
    res_sl <- phylolm(log.Sec1 ~ Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_sl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_sl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight.Mode)[1]
    coef = rbind(coef, coef_i)
  }
  return(coef)
}



## PART 2: loop through the 100 trees and merge all results
elev_bin_abs_100coef <- function(input, tag){
  
  coef_ntree = c()
  for (ii in 1:100){
    # use the new formulas
    coef_treeii = elev_bin_abs(treei = n_trees[[ii]], dataset = input)
    coef_treeii$treeID = ii
    coef_ntree = rbind(coef_ntree, coef_treeii)
  }
  saveRDS(coef_ntree, paste0('X:/home/_Wings_/HPC/elev_band_coef100/', tag, '.rds'))
  return(coef_ntree)
}


## results
elev_bin_abs_100coef(input = dat, tag = 'absolutue wing changes')



