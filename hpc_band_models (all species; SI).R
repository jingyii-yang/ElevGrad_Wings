setwd('/rdsgpfs/general/user/jy419/home/_Wings_/')

library(tidyverse)
{
  library(ape)
  library(caper)
  library(phylolm)
  library(mice)
}

n_trees <- readRDS('100 random trees.rds')
dat = read.csv('Analy data main.csv')


## make all categorical variables binary to compare effect sizes
{
  dat$Habitat.Openness = ifelse(dat$Habitat.Openness == 3, '1.open', '0.closed')
  dat$Migration = ifelse(dat$Migration == 3, '1.migratory', '0.sedentary')
  dat$AL.index = ifelse(dat$AL.index == 3, '1.frequent', '0.infrequent')
  
  dat$Flight_mode = ifelse(dat$Flight_mode == 'hover', 'flap', dat$Flight_mode)
  dat$Trophic.Level = ifelse(dat$Trophic.Level %in% c('Herbivore','Omnivore'), '1ry consumer', '2ry consumer')
}


dat$Elev = scale(dat$MAX)/2


## assign elevation bands to each species based on their max elevation
{
  dat$elevbin1 = ifelse(dat$MAX <= 3000, 1, 0) # the lowest MAX elev = 5m
  dat$elevbin2 = ifelse(dat$MAX >= 1000 & dat$MAX <= 4000, 1, 0)
  dat$elevbin3 = ifelse(dat$MAX >= 2000 & dat$MAX <= 5000, 1, 0)
  dat$elevbin4 = ifelse(dat$MAX >= 3000 & dat$MAX <= 6000, 1, 0)
  dat$elevbin5 = ifelse(dat$MAX >= 4000 & dat$MAX <= 7000, 1, 0)
  dat$elevbin6 = ifelse(dat$MAX >= 5000, 1, 0)
  # dat$elevbin6 = ifelse(dat$MAX >= 5000 & dat$MAX <= 8000, 1, 0) # only 1 species is above 8km (8300)
  # no difference including or excluding it.
}



## PART 1: do models within each elevation band, use the same tree for each of the 4 metrics

elev_bin <- function(treei, dataset){
  
  # choose the 6 elevation band columns
  bin_cols <<- str_which(colnames(dataset), 'elevbin')
  
  coef = c()
  for (i in bin_cols){
    
    ## select species in the current band (i)
    dat_bin = dataset %>% filter_at(i, all_vars(. == 1))
    
    datree_bin <- caper::comparative.data(data = dat_bin, phy = treei,          
                                          names.col = 'Species3',
                                          na.omit = FALSE)
    
    
    res_hwi <- phylolm(HWI ~ Elev * Flight_mode + Migration + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                       data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_hwi)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_hwi$formula[[2]])
    ## sample size depends on flight mode because adding the interaction means 
    ## the slope estimates represent the reference group only. 
    coef_i$sampleSize = table(dat_bin$Flight_mode)[1]  
    coef = rbind(coef, coef_i)
    
    
    res_wa <- phylolm(log.WA ~ Elev * Flight_mode + Migration + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wa)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wa$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight_mode)[1]
    coef = rbind(coef, coef_i)
    
    
    res_wl <- phylolm(log.Wing.Length ~ Elev * Flight_mode + Migration + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_wl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_wl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight_mode)[1]
    coef = rbind(coef, coef_i)
    
    
    res_sl <- phylolm(log.Sec1 ~ Elev * Flight_mode + Migration + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level,
                      data = datree_bin$data, phy = datree_bin$phy, model = 'lambda')
    coef_i = as.data.frame(summary(res_sl)$coefficients)
    coef_i$predictor = row.names(coef_i)
    coef_i$elevbin = i
    coef_i$Y = as.character(res_sl$formula[[2]])
    coef_i$sampleSize = table(dat_bin$Flight_mode)[1]
    coef = rbind(coef, coef_i)
    
  }
  return(coef)
}



## PART 2: loop through the 100 trees and merge all results

elev_bin_100coef <- function(input, tag){
  
  coef_ntree = c()
  for (ii in 1:length(n_trees)){
    coef_treeii = elev_bin(treei = n_trees[[ii]], dataset = input)
    coef_treeii$treeID = ii
    # this table will contain all coefficients from the 6 elevation bands x 4 wing metrics x 100 trees
    coef_ntree = rbind(coef_ntree, coef_treeii)
  }
    saveRDS(coef_ntree, paste0('/rdsgpfs/general/user/jy419/home/_Wings_/results/elev_band_coef100/', tag, '.rds'))
  return(coef_ntree)
}


########  MAIN ANALYSIS (GLOBAL MODEL) #########

# the full models (HWI, WA; WL, SL)
elev_bin_100coef(input = dat, tag = 'full')



######  Additional supplementary analysis 1: repeat with higher certainty data  ######

dat_A = filter(dat, AL_uncertainty == 'A') 

elev_bin_100coef(input = dat_A, tag = 'full high.certainty')




######  Additional supplementary analysis 2: repeat with mean elevation data  ######

dat$Elev = scale(dat$MEAN)/2 

### Note that the current elevation band division is still based on max elevation,
## need to re-divide the bands based on mean elevation to match the new model type.

{
  dat$elevbin1 = ifelse(dat$MEAN <= 3000, 1, 0) # the lowest MEAN elev = 2.5 m
  dat$elevbin2 = ifelse(dat$MEAN >= 1000 & dat$MEAN <= 4000, 1, 0)
  dat$elevbin3 = ifelse(dat$MEAN >= 2000 & dat$MEAN <= 5000, 1, 0)
  dat$elevbin4 = ifelse(dat$MEAN >= 3000, 1, 0)  # no species has a mean elevation > 6km (highest MEAN = 5100 m)
}

elev_bin_100coef(input = dat, tag = 'full mean.elevation')



