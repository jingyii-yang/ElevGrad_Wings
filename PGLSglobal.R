setwd("X:/home/_Wings_/HPC/")

library(tidyverse)
library(ape)
library(caper)
library(phylolm)
library(rr2)


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
  
  dat$TempVar = scale(log(dat$Temp..seasonality))/2
  dat$Lat = scale(sqrt(abs(dat$Centroid.latitude)))/2
}



########################################### MODELLING #####################################################

# Run each model 100 times to account for phylogenetic uncertainty, each time using a different tree.

res_ntree <- function(ntree=100, dataset, Y, mod_type, label){
  
  ## PART 1: run 100 models
  
  # raw model results to be saved in this list
  res = list()
  
  for (i in 1:ntree){
    
    datree_i <- caper::comparative.data(data = dataset, phy = n_trees[[i]],          
                                        names.col = 'Species_BirdTree',
                                        na.omit = FALSE)
    
    # all full models have the same set of predictors, only response variables (Y) are different
    if (mod_type==0) {
      res_i <- phylolm(formula(paste(Y, '~', 'Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level + log.Mass')),
                       data = datree_i$data, phy = datree_i$phy, model = 'lambda')}
    
    # although some supplementary analyses exclude body mass from predictors. Remove it.
    if (mod_type=='absolute_size') {
      res_i <- phylolm(formula(paste(Y, '~', 'Elev * Flight.Mode + Migration + Aerial.Lifestyle + Lat + TempVar + Habitat.Openness + Trophic.Level')),
                       data = datree_i$data, phy = datree_i$phy, model = 'lambda')}
    
    res[[i]] = res_i
  }
  
  # save raw model results for use later
  # saveRDS(res, paste0('../results/raw_models_100_global/', Y, label, '.rds'))
   saveRDS(res, paste0('X:/home/_Wings_/HPC/raw_models_100_global/', Y, label, '.rds'))
   
  
  
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
  
  saveRDS(coef, paste0('X:/home/_Wings_/HPC/global_coefs100/', Y, label, '.rds'))
  
  
  return(coef)
  
}



################################## Changes of relative wing metrics #############################

# Use maximum elevation as the main predictor
dat$Elev = scale(sqrt(dat$MAX))/2



######  Main analysis: all-species model [Figure3; Table S1 and S3 (1st rows)]
res_ntree(dataset = dat, Y = 'HWI',            label = '__HWI. all species', mod_type = 0) 
res_ntree(dataset = dat, Y = 'log.WA',          label = '__WA. all species', mod_type = 0)  
res_ntree(dataset = dat, Y = 'log.Wing.Length', label = '__WL. all species', mod_type = 0) 
res_ntree(dataset = dat, Y = 'log.Sec1',        label = '__SL. all species', mod_type = 0) 



######  Main supplementary analysis: non-migratory landbird model [Figure3; Table S1 and S3 (2nd rows)]
dat2 <- dat %>% filter(Migration != 3) %>% filter(Habitat.type != 'Marine')

res_ntree(dataset = dat2, Y = 'HWI',            label = '__HWI. landbird model', mod_type = 0) 
res_ntree(dataset = dat2, Y = 'log.WA',          label = '__WA. landbird model', mod_type = 0)  
res_ntree(dataset = dat2, Y = 'log.Wing.Length', label = '__WL. landbird model', mod_type = 0) 
res_ntree(dataset = dat2, Y = 'log.Sec1',        label = '__SL. landbird model', mod_type = 0) 



      ##### (set soaring species as the reference level to get direct estimates for soaring species) 
      #     for the extra predictor line in Table S1 and S3 >> Elevation (soaring) *
      dat_soar = dat
      dat_soar$Flight.Mode[dat_soar$Flight.Mode=='soaring'] <- '0.soaring'
      dat_soar2 <- dat_soar %>% filter(Migration != 3) %>% filter(Habitat.type != 'Marine')
      
      # All-species model
      res_ntree(dataset = dat_soar, Y = 'HWI',            label = '__HWI. all species (ref=soaring)', mod_type = 0)  
      res_ntree(dataset = dat_soar, Y = 'log.WA',          label = '__WA. all species (ref=soaring)', mod_type = 0)  
      res_ntree(dataset = dat_soar, Y = 'log.Wing.Length', label = '__WL. all species (ref=soaring)', mod_type = 0) 
      res_ntree(dataset = dat_soar, Y = 'log.Sec1',        label = '__SL. all species (ref=soaring)', mod_type = 0) 
      
      # Non-migratory landbird model
      res_ntree(dataset = dat_soar2, Y = 'HWI',            label = '__HWI. landbird model (ref=soaring)', mod_type = 0) 
      res_ntree(dataset = dat_soar2, Y = 'log.WA',          label = '__WA. landbird model (ref=soaring)', mod_type = 0)  
      res_ntree(dataset = dat_soar2, Y = 'log.Wing.Length', label = '__WL. landbird model (ref=soaring)', mod_type = 0) 
      res_ntree(dataset = dat_soar2, Y = 'log.Sec1',        label = '__SL. landbird model (ref=soaring)', mod_type = 0) 



######  Additional supplementary analysis 1: repeat with higher certainty data [Table S1 and S3 (3rd rows)]

dat_A = filter(dat, Aerial.Lifestyle.Certainty == 'A') 

res_ntree(dataset = dat_A, Y = 'HWI',            label = '__HWI high certainty', mod_type = 0)  
res_ntree(dataset = dat_A, Y = 'log.WA',          label = '__WA high certainty', mod_type = 0)  
res_ntree(dataset = dat_A, Y = 'log.Wing.Length', label = '__WL high certainty', mod_type = 0) 
res_ntree(dataset = dat_A, Y = 'log.Sec1',        label = '__SL high certainty', mod_type = 0) 



######  Additional supplementary analysis 2: repeat with mean elevation data [Table S1 and S3 (4th rows)]

dat_mean = dat
dat_mean$Elev = scale(sqrt(dat_mean$MEAN))/2

res_ntree(dataset = dat_mean, Y = 'HWI',            label = '__HWI mean elevation', mod_type = 0)  
res_ntree(dataset = dat_mean, Y = 'log.WA',          label = '__WA mean elevation', mod_type = 0)  
res_ntree(dataset = dat_mean, Y = 'log.Wing.Length', label = '__WL mean elevation', mod_type = 0) 
res_ntree(dataset = dat_mean, Y = 'log.Sec1',        label = '__SL mean elevation', mod_type = 0) 



#################################  Changes of body mass #########################################


# All-species model [Table S2 (1st row)]
res_ntree(dataset = dat, Y = 'log.Mass', label = '__abs_Mass. full', mod_type = 'absolute_size') 


# Non-migratory landbird model [Table S2 (2nd row)]
res_ntree(dataset = dat2, Y = 'log.Mass', label = '__abs_Mass. landbird', mod_type = 'absolute_size')  



################################## Changes of absolute wing metrics #############################


# All-species model [Table S1 and S3 (5th rows)]
res_ntree(dataset = dat, Y = 'HWI',            label = '__abs_HWI full', mod_type = 'absolute_size') 
res_ntree(dataset = dat, Y = 'log.WA',          label = '__abs_WA full', mod_type = 'absolute_size')  
res_ntree(dataset = dat, Y = 'log.Wing.Length', label = '__abs_WL full', mod_type = 'absolute_size') 
res_ntree(dataset = dat, Y = 'log.Sec1',        label = '__abs_SL full', mod_type = 'absolute_size') 


      ##### (set soaring species as the reference level to get direct estimates for soaring species) 
      #     for the extra predictor line in Table S1 and S3 >> Elevation (soaring) *
      res_ntree(dataset = dat_soar, Y = 'HWI',            label = '__abs_HWI. all species (ref=soaring)', mod_type = 'absolute_size')  
      res_ntree(dataset = dat_soar, Y = 'log.WA',          label = '__abs_WA. all species (ref=soaring)', mod_type = 'absolute_size')  
      res_ntree(dataset = dat_soar, Y = 'log.Wing.Length', label = '__abs_WL. all species (ref=soaring)', mod_type = 'absolute_size') 
      res_ntree(dataset = dat_soar, Y = 'log.Sec1',        label = '__abs_SL. all species (ref=soaring)', mod_type = 'absolute_size') 
      
