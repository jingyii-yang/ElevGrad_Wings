library(tidyverse)
library(ggtext)


######################## Figure S2. correlation between HWA and total aeronautical wing area ####################

# merge cleaned literature data on total wing area

library(tidyverse)
setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/source data (elev, wing area)/cleaned wing area data/")
{
indiv = read.csv('indiv sp.csv')
bru_data = read.csv('cleaned Bruderer et al. 2010.csv')
clm1_data = read.csv('cleaned Claramunt 2012.csv')
clm2_data = read.csv('cleaned Claramunt 2021.uk.csv')
hwd_data = read.csv('cleaned howard data 2018.csv')
malo_data = read.csv('cleaned malo_data_2021.csv')
wtn_data = read.csv('cleaned Watanabe data 2016.csv')
fu1 = read.csv('cleaned Fu 2023 (1).csv')
fu2 = read.csv('cleaned Fu 2023 (2).csv')
srn = read.csv('cleaned Serrano 2016.csv')
vsc = read.csv('cleaned Viscor 1987.csv')

# add reference info
bru_data$ref = 'Bruderer et al. 2010'
clm1_data$ref = 'Claramunt et al. 2012'
clm2_data$ref = 'Claramunt et al. 2021'
hwd_data$ref = 'Howard et al. 2018'
malo_data$ref = 'Malo et al. 2021'
wtn_data$ref = 'Watanabe 2016'
fu1$ref = 'Fu et al. 2023 (A1)'
fu2$ref = 'Fu et al. 2023 (A2)'
srn$ref = 'Serrano et al. 2016'
vsc$ref = 'Viscor et al. 1987'
}

# standardize units: wing area to mm2
{
bru_data$wing.area..m2.= bru_data$wing.area..m2. * 1000000
  
clm1_data$Total.wing.area.cm2 = clm1_data$Total.wing.area.cm2 * 100

clm2_data$wing.area.m2 = clm2_data$wing.area.m2 * 1000000

hwd_data$Wing.area.m2 = hwd_data$Wing.area.m2 * 1000000

malo_data$Wing.area..S..in.m2.= malo_data$Wing.area..S..in.m2.* 1000000

wtn_data$Wing.area.m2 = wtn_data$Wing.area.m2 * 1000000

fu1$TotalWingArea.cm2 = fu1$TotalWingArea.cm2 * 100

fu2$Wing.area.calculated.cm2 = fu2$Wing.area.calculated.cm2 * 100

srn$log.SL.lift.surface.m2 = 10^(srn$log.SL.lift.surface.m2) * 1000000  # raw data on log scale

vsc$lift.wing.area.cm2 = vsc$lift.wing.area.cm2 * 100
}

# standardise column names and merge datasets
{
coln = c('Species', 'Wing.span.mm', 'Wing.area.mm2')
colnames(indiv)[c(1,2,3)] <- coln
colnames(bru_data)[c(1,2,3)] <- coln
colnames(clm1_data)[c(1,2,3)] <- coln
colnames(clm2_data)[c(1,3,4)] <- coln
colnames(hwd_data)[c(1,2,3)] <- coln
colnames(malo_data)[c(1,5,3)] <- coln
colnames(wtn_data)[c(1,2,3)] <- coln
colnames(fu1)[c(1,2,3)] <- coln
colnames(fu2)[c(1,2,3)] <- coln
colnames(srn)[c(1,2,3)] <- coln
colnames(vsc)[c(1,3,4)] <- coln
}

coln = c('Species', 'Wing.span.mm', 'Wing.area.mm2', 'ref')
Wing_Data = rbind(indiv[,coln], bru_data[,coln],
                  clm1_data[,coln],clm2_data[,coln],
                  hwd_data[,coln], malo_data[,coln], wtn_data[,coln],
                  fu1[,coln], fu2[,coln], srn[,coln], vsc[,coln])


# clean species names before calculating the mean
{
# remove special characters
Wing_Data$Species = str_trim(Wing_Data$Species, side = 'both')
Wing_Data$Species = janitor::make_clean_names(Wing_Data$Species, case = 'sentence')
Wing_Data$Species = gsub(' ', '_', Wing_Data$Species)

# remove subspecies info
Wing_Data$Genus <- sapply(str_split(Wing_Data$Species, '_'), '[[', 1)
Wing_Data$sp <- sapply(str_split(Wing_Data$Species, '_'), '[[', 2)
Wing_Data$Species <- paste(Wing_Data$Genus, Wing_Data$sp, sep = '_')

## align species names with BirdTree names
sp_match = read.csv('sp names mismatches.csv')
Wing_Data = left_join(Wing_Data, sp_match, by=c('Species' = 'names'))
Wing_Data$jetz.name = ifelse(is.na(Wing_Data$jetz.name), Wing_Data$Species, Wing_Data$jetz.name)
}


## create total wing area dataset
Wing_Data = Wing_Data %>% group_by(jetz.name) %>% summarise(Wing.area.mm2 = mean(Wing.area.mm2, na.rm = T),
                                                            ref = unique(ref)[1])

dat <- read.csv("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/data/Analy data main.csv")
test=left_join(dat, Wing_Data, by=c('Species3'='jetz.name'))
  colSums(!is.na(test))



# calculate HWA (both hand portions) using formula by Wright et al. 2014
test$WA.hand = test$Wing.Length*test$Secondary1*pi/2

test_area = na.omit(test[,c('Wing.area.mm2','WA.hand')]) %>% log()
cor1 = cor.test(test_area$Wing.area.mm2, test_area$WA.hand)
    
pdf('../../figures/SI/Fig. S2 correlation between total and hand wing area (log, scatter).pdf', width = 5, height = 5)
ggplot(test_area, aes(Wing.area.mm2, WA.hand))+theme_classic()+
    geom_point(size = 1, alpha = 0.5) +
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=14)) +
    labs(x='Total wing area (log)', y = 'Hand-wing area (log)') +
    annotate('text', label=paste0("Pearson's correlation = ", round(cor1$estimate, 2),
                                  '\np < 0.001 \nn = ', nrow(test_area), ' species'), 
             hjust=0,
             x = 7, y = 12.5, size = 4.7)

dev.off()




########### Fig S7-8: correlation between different elevation data methods ###############

setwd('C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/data/')
elev = read.csv('RAW elev data year around data.csv')


####..... Fig.S7. correlation between raw data collected from different ranges ############

# merge data from ref. del Hoyo (2020) and ref. White et al. (2015), both extracted data based on all known ranges of the species
elev$all_ranges_max = ifelse(is.na(elev$new_max), elev$WB_max, elev$new_max)
elev$all_ranges_min = ifelse(is.na(elev$new_min), elev$WB_min, elev$new_min)
# add mean elevation
elev$all_ranges_mean = (elev$all_ranges_max + elev$all_ranges_min)/2
elev$QJ_mean.breeding. = (elev$QJ_max.breeding. + elev$QJ_min.breeding.)/2


pdf('../../elev-r/figures/SI/Fig. S7 Correlation between raw elevation data (breeding vs combined all ranges; scatter).pdf', width = 11, height = 5)

elev_plot = na.omit(elev[, c('all_ranges_max', 'QJ_max.breeding.', 'Migration')])
elev_plot$Migration = ifelse(elev_plot$Migration == 3, 'Yes', 'No')
cor2 = cor.test(elev_plot$all_ranges_max, elev_plot$QJ_max.breeding.)
pa=ggplot(elev_plot, aes(all_ranges_max/1000, QJ_max.breeding./1000, col=Migration))+theme_classic()+
          geom_point(size = 1, alpha = 0.5) +
          theme(axis.title = element_text(size=16),
                axis.text = element_text(size=16), legend.position = 'none') +
          labs(x='Elevation - all ranges (km)', y = 'Elevation - breeding range (km)') +
          scale_x_continuous(limits = c(0, 8.500)) +
      scale_color_manual(values = c('Black', 'red')) +
          annotate('text', label=paste0("Pearson's correlation = ", round(cor2$estimate, 2),
                                '\np < 0.001 \nn = ', nrow(elev_plot), ' species'), 
           hjust=0, x = 0, y = 7.700, size = 4.8)

elev_plot = na.omit(elev[, c('all_ranges_mean', 'QJ_mean.breeding.', 'Migration')])
elev_plot$Migration = ifelse(elev_plot$Migration == 3, 'Yes', 'No')
cor3 = cor.test(elev_plot$all_ranges_mean, elev_plot$QJ_mean.breeding.)
pb=ggplot(elev_plot, aes(all_ranges_mean/1000, QJ_mean.breeding./1000, col=Migration))+theme_classic()+
          geom_point(size = 1, alpha = 0.5) +
          theme(axis.title = element_text(size=16),
                axis.text = element_text(size=16), legend.position = 'none') +
          labs(x='Elevation - all ranges (km)', y = 'Elevation - breeding range (km)') +
          scale_x_continuous(limits = c(0, 6.200)) +
          scale_y_continuous(limits = c(0, 6.200)) +
      scale_color_manual(values = c('Black', 'red')) +
          annotate('text', label=paste0("Pearson's correlation = ", round(cor3$estimate, 2),
                                        '\np < 0.001 \nn = ', nrow(elev_plot), ' species'), 
           hjust=0, x = 0, y = 5.800, size = 4.8)

cowplot::plot_grid(pa, pb, scale = 0.97)
dev.off()


## info on the % of non-migratory/migratory species
table(elev$Migration[!is.na(elev$all_ranges_max) & !is.na(elev$QJ_max.breeding.)])
table(elev$Migration[!is.na(elev$all_ranges_mean) & !is.na(elev$QJ_mean.breeding.)])



####..... Fig.S8 correlation between final elevation data selected using methods 1 & 2 ####

## select elevation data: 
# Method 1: prioritise newer data sources that extracted data based on widest ranges
{
  elev$MAX <- elev$new_max
  elev$MAX <- ifelse(is.na(elev$MAX), elev$WB_max, elev$MAX)
  elev$MAX <- ifelse(is.na(elev$MAX), elev$QJ_max.breeding., elev$MAX)
  
  elev$MIN <- elev$new_min
  elev$MIN <- ifelse(is.na(elev$MIN), elev$WB_min, elev$MIN)
  elev$MIN <- ifelse(is.na(elev$MIN), elev$QJ_min.breeding., elev$MIN)
  
  elev$MEAN = (elev$MAX + elev$MIN)/2
}

# Method 2: use the extreme values of the max and min elevation reported for species, with no preference to sources.
{
  elev$MAX1 <- apply(elev[, c('new_max', 'WB_max', 'QJ_max.breeding.')], 1, max, na.rm=T)
  elev$MIN1 <- apply(elev[, c('new_min', 'WB_min', 'QJ_min.breeding.')], 1, min, na.rm=T)
  elev$MEAN1 = (elev$MAX1 + elev$MIN1)/2
  
  elev$MAX1 <- as.numeric(elev$MAX1)
  elev$MIN1 <- as.numeric(elev$MIN1)
  elev$MEAN1 <- as.numeric(elev$MEAN1)
}


pdf('../../elev-r/figures/SI/Fig.S8 correlation between elevation measures (max 1-2, mean 1-2).pdf', width = 11, height = 5)

elev_plot = na.omit(elev[, c('MAX', 'MAX1')])
cor4 = cor.test(elev_plot$MAX, elev_plot$MAX1)
pa=ggplot(elev_plot, aes(MAX/1000, MAX1/1000))+theme_classic()+
          geom_point(size = 1, alpha = 0.5) +
          theme(axis.title = element_text(size=16),
                axis.text = element_text(size=16)) +
          labs(x='Maximum elevation 1 (km)', y = 'Maximum elevation 2 (km)') +
          scale_x_continuous(limits = c(0, 8.500)) +
          annotate('text', label=paste0("Pearson's correlation = ", round(cor4$estimate, 2),
                                        '\np < 0.001 \nn = 9986 species'), 
                   hjust=0, x = 0, y = 7.700, size = 4.8)

elev_plot = na.omit(elev[, c('MEAN', 'MEAN1')])
cor5 = cor.test(elev_plot$MEAN, elev_plot$MEAN1)
pb=ggplot(elev_plot, aes(MEAN/1000, MEAN1/1000))+theme_classic()+
          geom_point(size = 1, alpha = 0.5) +
          theme(axis.title = element_text(size=16),
                axis.text = element_text(size=16)) +
          labs(x='Mean elevation 1 (km)', y = 'Mean elevation 2 (km)') +
          scale_x_continuous(limits = c(0, 6.200)) +
          scale_y_continuous(limits = c(0, 6.200)) +
          annotate('text', label=paste0("Pearson's correlation = ", round(cor5$estimate, 2),
                                        '\np < 0.001 \nn = 9986 species'), 
                   hjust=0, x = 0, y = 5.800, size = 4.8)

cowplot::plot_grid(pa, pb, scale = 0.97)
dev.off()

