setwd('/rdsgpfs/general/user/jy419/home/_Wings_/')

library(tidyverse)
{
  library(ape)
  library(caper)
  library(phylolm)
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




res_ntree <- function(ntree=100, dataset, Y, mod_type, label){
  
  ## PART 1: run 100 models
  
  # raw model results to be saved in this list
  res = list()
  # run each model 100 times to account for phylogenetic uncertainty
  for (i in 1:ntree){
    datree_i <- caper::comparative.data(data = dataset, phy = n_trees[[i]],          
                                        names.col = 'Species3',
                                        na.omit = FALSE)
    
    # all full models have the same set of predictors, only response variables (Y) are different
    if (mod_type==0) {
      res_i <- phylolm(formula(paste(Y, '~', 'Elev * Flight_mode + Migration + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level')),
                       data = datree_i$data, phy = datree_i$phy, model = 'lambda')}
    
    # although some supplementary analyses require reduced set of predictors. Remove those predictors.
    if (mod_type=='landbird') {
      res_i <- phylolm(formula(paste(Y, '~', 'Elev * Flight_mode + AL.index + Lat + TempVar + Habitat.Openness + log.Mass + Trophic.Level')),
                       data = datree_i$data, phy = datree_i$phy, model = 'lambda')}
    
    res[[i]] = res_i
  }
  
  
  saveRDS(res, paste0('/rdsgpfs/general/user/jy419/home/_Wings_/results/raw_model_100_global/', Y, label, '.rds'))
  
  
  
  ## PART 2: merge all results
  
  coef=c()
  # format results tables so they can be used in the pool function later
  for (z in 1:length(res) ){
    coef_i <- as.data.frame(summary(res[[z]])$coefficients[,c(1:2)]) 
    coef_i$term = row.names(coef_i)
    coef_i$df.residual = length(res[[z]]$y) - length(res[[z]]$coefficients)
    coef = rbind(coef, coef_i)
  }
  colnames(coef)[c(1,2)] <- c('estimate', 'std.error')
  
  saveRDS(coef, paste0('/rdsgpfs/general/user/jy419/home/_Wings_/results/global_coefs100/', Y, label, '.rds'))
  
  
  return(coef)
  
}



########  MAIN ANALYSIS (GLOBAL MODEL) #########

### each model takes ~40 mins to run on an average PC.

# All-species model
res_ntree(dataset = dat, Y = 'HWI', label = '__HWI (ref=flap) full', mod_type = 0) 
res_ntree(dataset = dat, Y = 'log.WA', label = '__WA (ref=flap) full', mod_type = 0)  
res_ntree(dataset = dat, Y = 'log.Wing.Length', label = '__WL (ref=flap) full', mod_type = 0) 
res_ntree(dataset = dat, Y = 'log.Sec1', label = '__SL (ref=flap) full', mod_type = 0) 


# Main sensitivity analysis: Non-migratory landbird model
dat2 <- dat %>% filter(Migration == '0.sedentary') %>% filter(Habitat != 'Marine')

res_ntree(dataset = dat2, Y = 'HWI', label = '__HWI (ref=flap) land.sedentary', mod_type = 'landbird') 
res_ntree(dataset = dat2, Y = 'log.WA', label = '__WA (ref=flap) land.sedentary', mod_type = 'landbird')  
res_ntree(dataset = dat2, Y = 'log.Wing.Length', label = '__WL (ref=flap) land.sedentary', mod_type = 'landbird') 
res_ntree(dataset = dat2, Y = 'log.Sec1', label = '__SL (ref=flap) land.sedentary', mod_type = 'landbird') 



      ##### (set soaring species as the reference level to get direct estimates for soaring species) #####
      
      dat$Flight_mode[dat$Flight_mode=='soar'] <- '0.soar'
      
      # All-species model
      res_ntree(dataset = dat, Y = 'HWI', label = '__HWI (ref=soar) full', mod_type = 0)  
      res_ntree(dataset = dat, Y = 'log.WA', label = '__WA (ref=soar) full', mod_type = 0)  
      res_ntree(dataset = dat, Y = 'log.Wing.Length', label = '__WL (ref=soar) full', mod_type = 0) 
      res_ntree(dataset = dat, Y = 'log.Sec1', label = '__SL (ref=soar) full', mod_type = 0) 
      
      # Main sensitivity analysis: Non-migratory landbird model
      dat2 <- dat %>% filter(Migration == '0.sedentary') %>% filter(Habitat != 'Marine')
      
      res_ntree(dataset = dat2, Y = 'HWI', label = '__HWI (ref=soar) land.sedentary', mod_type = 'landbird') 
      res_ntree(dataset = dat2, Y = 'log.WA', label = '__WA (ref=soar) land.sedentary', mod_type = 'landbird')  
      res_ntree(dataset = dat2, Y = 'log.Wing.Length', label = '__WL (ref=soar) land.sedentary', mod_type = 'landbird') 
      res_ntree(dataset = dat2, Y = 'log.Sec1', label = '__SL (ref=soar) land.sedentary', mod_type = 'landbird') 





######  Additional supplementary analysis 1: repeat with higher certainty data  ######

dat_A = filter(dat, AL_uncertainty == 'A') 

res_ntree(dataset = dat_A, Y = 'HWI', label = '__HWI (ref=flap) full.high certainty', mod_type = 0)  
res_ntree(dataset = dat_A, Y = 'log.WA', label = '__WA (ref=flap) full.high certainty', mod_type = 0)  
res_ntree(dataset = dat_A, Y = 'log.Wing.Length', label = '__WL (ref=flap) full.high certainty', mod_type = 0) 
res_ntree(dataset = dat_A, Y = 'log.Sec1', label = '__SL (ref=flap) full.high certainty', mod_type = 0) 



######  Additional supplementary analysis 2: repeat with mean elevation data  ######

dat_mean = dat
dat_mean$Elev = scale(dat_mean$MEAN)/2

res_ntree(dataset = dat_mean, Y = 'HWI', label = '__HWI (ref=flap) full.mean elevation', mod_type = 0)  
res_ntree(dataset = dat_mean, Y = 'log.WA', label = '__WA (ref=flap) full.mean elevation', mod_type = 0)  
res_ntree(dataset = dat_mean, Y = 'log.Wing.Length', label = '__WL (ref=flap) full.mean elevation', mod_type = 0) 
res_ntree(dataset = dat_mean, Y = 'log.Sec1', label = '__SL (ref=flap) full.mean elevation', mod_type = 0) 

