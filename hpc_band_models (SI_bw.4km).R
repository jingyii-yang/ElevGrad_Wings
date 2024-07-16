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
  dat$Trophic.Level = ifelse(dat$Trophic.Level %in% c('Herbivore','Omnivore'), '1ry consumer', '2ry consumer')
}


dat$Elev = scale(dat$MAX)/2


## assign elevation bands to each species
{
  dat$elevbin1 = ifelse(dat$MAX <= 4000, 1, 0) # the lowest MAX elev = 5m
  dat$elevbin2 = ifelse(dat$MAX >= 1000 & dat$MAX <= 5000, 1, 0)
  dat$elevbin3 = ifelse(dat$MAX >= 2000 & dat$MAX <= 6000, 1, 0)
  dat$elevbin4 = ifelse(dat$MAX >= 3000 & dat$MAX <= 7000, 1, 0)
  dat$elevbin5 = ifelse(dat$MAX >= 4000, 1, 0)
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
    ## the slope estimate can represent for the reference group only. 
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
    # this table will contain all coefficients from the 5 elevation bands x 4 wing metrics x 100 trees
    coef_ntree = rbind(coef_ntree, coef_treeii)
  }
    saveRDS(coef_ntree, paste0('/rdsgpfs/general/user/jy419/home/_Wings_/results/SI4_elev_band_coef100_BW.4km/', tag, '.rds'))
  return(coef_ntree)
}


########  MAIN ANALYSIS (GLOBAL MODEL) #########

# the full models (HWI, WA; WL, SL)
elev_bin_100coef(input = dat, tag = 'full (bw=4km)')

