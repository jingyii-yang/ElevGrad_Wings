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


## make categorical variables binary to compare effect sizes
{
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



######  Supplementary analysis 1: non-migratory landbirds  ###############

dat2 <- dat %>% filter(Migration != 3) %>% filter(Habitat != 'Marine')
elev_bin_100coef(input = dat2, tag = 'land.sedentary')



######  Supplementary analysis 2: retest with high certainty data  ######

dat_A = filter(dat, AL_uncertainty == 'A') 
elev_bin_100coef(input = dat_A, tag = 'full high.certainty')



######  Supplementary analysis 3: retest with mean elevation data  ######

dat_mean = dat
dat_mean$Elev = scale(sqrt(dat_mean$MEAN))/2 

## need to re-divide the bands based on mean elevation to match the new model structure.
{# Remove existing 'elevbin' columns first
  dat_mean = dplyr::select(dat_mean, -contains('elevbin'))
  
  dat$elevbin1 = ifelse(dat$MEAN <= 3000, 1, 0) # the lowest MEAN elev = 2.5 m
  dat$elevbin2 = ifelse(dat$MEAN >= 1000 & dat$MEAN <= 4000, 1, 0)
  dat$elevbin3 = ifelse(dat$MEAN >= 2000 & dat$MEAN <= 5000, 1, 0)
  dat$elevbin4 = ifelse(dat$MEAN >= 3000, 1, 0)  # no species has a mean elevation > 6km (highest MEAN = 5100 m)
}

elev_bin_100coef(input = dat, tag = 'full mean.elevation')



######  Supplementary analysis 4.1: use bandwidth = 2km  ##############

dat_2km = dat

## re-assign elevation bands to each species

{ # Remove existing 'elevbin' columns
  dat_2km = dplyr::select(dat_2km, -contains('elevbin'))
  
  dat_2km$elevbin1 = ifelse(dat_2km$MAX <= 2000, 1, 0) # the lowest MAX elev = 5m
  dat_2km$elevbin2 = ifelse(dat_2km$MAX >= 1000 & dat_2km$MAX <= 3000, 1, 0)
  dat_2km$elevbin3 = ifelse(dat_2km$MAX >= 2000 & dat_2km$MAX <= 4000, 1, 0)
  dat_2km$elevbin4 = ifelse(dat_2km$MAX >= 3000 & dat_2km$MAX <= 5000, 1, 0)
  dat_2km$elevbin5 = ifelse(dat_2km$MAX >= 4000 & dat_2km$MAX <= 6000, 1, 0)
  dat_2km$elevbin6 = ifelse(dat_2km$MAX >= 5000 & dat_2km$MAX <= 7000, 1, 0)
  #  dat_2km$elevbin7 = ifelse(dat_2km$MAX >= 6000, 1, 0) this band contains too few species (18) so the model will not work.
}
elev_bin_100coef(input = dat_2km, tag = 'full (bw=2km)')



######  Supplementary analysis 4.2: use bandwidth = 4km  ###############

dat_4km = dat

## re-assign elevation bands to each species

{ # Remove existing 'elevbin' columns
  dat_4km = dplyr::select(dat_4km, -contains('elevbin'))
  
  dat_4km$elevbin1 = ifelse(dat_4km$MAX <= 4000, 1, 0) # the lowest MAX elev = 5m
  dat_4km$elevbin2 = ifelse(dat_4km$MAX >= 1000 & dat_4km$MAX <= 5000, 1, 0)
  dat_4km$elevbin3 = ifelse(dat_4km$MAX >= 2000 & dat_4km$MAX <= 6000, 1, 0)
  dat_4km$elevbin4 = ifelse(dat_4km$MAX >= 3000 & dat_4km$MAX <= 7000, 1, 0)
  dat_4km$elevbin5 = ifelse(dat_4km$MAX >= 4000, 1, 0)
}
elev_bin_100coef(input = dat_4km, tag = 'full (bw=4km)')


