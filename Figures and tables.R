setwd("X:/home/_Wings_/HPC/")

{
library(tidyverse)
library(ggtext)
library(cowplot)

library(ape)
  library(caper)
library(phylolm)
  library(mice)
library(rr2)
}

select = dplyr::select
between = dplyr::between

# read the dataset for main analyses
dat = read.csv('Data S1_BirdTree.csv') %>% rename(MAX = 'Max.Elevation.1', MEAN = 'Mean.Elevation.1') %>%
  filter(Flight.Mode != 'flightless')


##### Data transformation (same as models)
{
  dat$log.Mass <- scale(log(dat$Body.Mass))/2
  
  dat$log.Wing.Length <- scale(log(dat$Wing.Length))/2
  dat$log.Sec1 <- scale(log(dat$Secondary1))/2
  dat$log.WA = scale(log(dat$Wing.Length * dat$Secondary1 * pi/2))/2
  dat$HWI <- scale(dat$Hand.Wing.Index)/2
  
  dat$TempVar = scale(log(dat$Temp..seasonality))/2
  dat$Lat = scale(sqrt(abs(dat$Centroid.latitude)))/2
}

## calculate relative wing metrics for plotting
dat$rela.HWI = lm(dat$HWI ~ dat$log.Mass)$residuals
dat$rela.WA = lm(dat$log.WA ~ dat$log.Mass)$residuals
    

############## ......................... Figure 2 .........................##############

######## top panels: scatter plots + density plots ########

cus_cols = c(# for the main scatter plot
             '#bedee8', '#012a7a', 
             # regression line would be invisible if using exactly same colour as the points,
             # use darker shades
             '#669eb0', '#153a85',
             # use matching SE shades
             '#669eb0', '#435d99',
             # density plots overlapping so they must use transparent and darker colours 
             '#53d3fc', '#020a61')
            

p1=ggplot(dat, aes(MAX, rela.HWI, col=Flight.Mode)) + geom_point(alpha=0.4) +
    theme_classic() +
    scale_colour_manual(values = c(cus_cols[1], cus_cols[2])) + 
    # set the x axis limits and expand to align with side panels later.
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000), labels = c(0, 2, 4, 6, 8), limits = c(-400, 8700), expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
    labs(x = 'Maximum elevation (km)', y = 'Relative hand-wing index', col = 'Flight mode') +
    # regression line of soaring species
    geom_smooth(data = dat[dat$Flight.Mode == 'soaring', ], col=cus_cols[4], fill=cus_cols[6], method = 'gam') +
    # remove marine soaring species (dotted line)
    geom_smooth(data = dat[dat$Flight.Mode == 'soaring' & dat$Habitat.type != 'Marine', ], col=cus_cols[4], se=F, method = 'gam', linetype = 2, linewidth = 0.7) +
    # regression line of flapping species
    geom_smooth(data = dat[dat$Flight.Mode == 'flapping', ], col=cus_cols[3], fill=cus_cols[5], method = 'gam') +
        theme(legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16))
           

p2=ggplot(dat, aes(MAX, rela.WA, col=Flight.Mode)) + geom_point(alpha=0.4) +
    theme_classic() +
    scale_colour_manual(values = c(cus_cols[1], cus_cols[2])) + 
    # set the x axis limits and expand to align with side panels later.
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000), labels = c(0, 2, 4, 6, 8), limits = c(-400, 8700), expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(x = 'Maximum elevation (km)', y = 'Relative hand-wing area', col = 'Flight mode') +
    # regression line of soaring species
    geom_smooth(data = dat[dat$Flight.Mode == 'soaring', ], col=cus_cols[4], fill=cus_cols[6], method = 'gam') +
    # remove marine soaring species (dotted line)
    geom_smooth(data = dat[dat$Flight.Mode == 'soaring' & dat$Habitat.type != 'Marine', ], col=cus_cols[4], se=F, method = 'gam', linetype = 2, linewidth = 0.7) +
    # regression line of flapping species
    geom_smooth(data = dat[dat$Flight.Mode == 'flapping', ], col=cus_cols[3], fill=cus_cols[5], method = 'gam') +
      # add legend
    annotate(geom = "rect", xmin = 6400, xmax = 7000, ymin = -0.56, ymax = -0.44,
             fill = '#bedee8', colour = "white", alpha = 0.9) + 
    annotate(geom = "rect", xmin = 6400, xmax = 7000, ymin = -0.68, ymax = -0.56,
             fill = '#153a85', colour = "white", alpha = 0.9) +  
    annotate(geom = "text", x = 7800, y = -0.5, label = 'Flapping') +
    annotate(geom = "text", x = 7720, y = -0.62, label = 'Soaring') +
    theme(legend.position = 'none',
            axis.title = element_text(size = 18),
              axis.text = element_text(size = 16))


# make side density plots for panel a.

xdens1 <- axis_canvas(p1, axis = "x") +
    geom_density(data = dat, aes(x = MAX, fill = Flight.Mode),
                alpha = 0.3, size = 0.25) +
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    # to align the x-axis of side panels to the main scatter plot.   
    scale_x_continuous(limits = c(-400, 8700))

ydens1 <- axis_canvas(p1, axis = "y", coord_flip = TRUE) +
    geom_density(data = dat, aes(x = rela.HWI, fill = Flight.Mode),
                  alpha = 0.3, size = 0.25)+
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    coord_flip()

# figure 2A
t1 <- insert_xaxis_grob(p1, xdens1, grid::unit(.1, "null"), position = "top")
t2<- insert_yaxis_grob(t1, ydens1, grid::unit(.1, "null"), position = "right")


# make side density plots for panel b.
xdens2 <- axis_canvas(p2, axis = "x") +
    geom_density(data = dat, aes(x = MAX, fill = Flight.Mode),
                alpha = 0.3, size = 0.25) +
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    # to align the x-axis of side panels to the main scatter plot.   
    scale_x_continuous(limits = c(-400, 8700))

ydens2 <- axis_canvas(p2, axis = "y", coord_flip = TRUE) +
    geom_density(data = dat, aes(x = rela.WA, fill = Flight.Mode),
                  alpha = 0.3, size = 0.25)+
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    coord_flip()

# figure 2B
t3 <- insert_xaxis_grob(p2, xdens2, grid::unit(.1, "null"), position = "top")
t4<- insert_yaxis_grob(t3, ydens2, grid::unit(.1, "null"), position = "right")



########  bottom panels: whisker plots  ########

datbin = dat
# divide species into 1km-wide elevational bands
datbin$elevbin = datbin$MAX %/% 1000 
# calculate data distribution of HWI and HWA for each band
datbin.plot <- datbin %>% group_by(elevbin) %>% summarise(mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                            lower.hwi = quantile(rela.HWI, 0.25, na.rm=TRUE),
                                                            upper.hwi = quantile(rela.HWI, 0.75, na.rm=TRUE),
                                                          mid.wa = median(rela.WA, na.rm = TRUE),
                                                            lower.wa = quantile(rela.WA, 0.25, na.rm=TRUE),
                                                            upper.wa = quantile(rela.WA, 0.75, na.rm=TRUE),
                                                          n=n())

# add the regression coefficients
stat.c = paste0('Slope = ',  round(lm(data = datbin.plot, mid.hwi ~ elevbin)$coefficients[2], 3), ', ', '*p*', ' < 0.001',
                          '<br> *R*<sup>2</sup> = ',  round(summary(lm(data = datbin.plot, mid.hwi ~ elevbin))$adj.r.squared, 2))

stat.d = paste0('Slope = ',  round(lm(data = datbin.plot, mid.wa ~ elevbin)$coefficients[2], 3), ', ', '*p*', ' < 0.01',
                          '<br> *R*<sup>2</sup> = ',  round(summary(lm(data = datbin.plot, mid.wa ~ elevbin))$adj.r.squared, 2))
  

q1= ggplot(datbin.plot, aes(elevbin, mid.hwi)) + 
  # median HWI 
    geom_point(aes(size = n), col = 'brown') +
  # 1st and 3rd quartiles of HWI
    geom_errorbar(aes(ymin = lower.hwi, ymax = upper.hwi), col = 'grey50') + 
      scale_y_continuous(limits = c(-0.7, 1), breaks = c(-0.5, 0, 0.5)) +
      scale_x_continuous(breaks = c(0,2,4,6,8), labels = c('0-1','2-3', '4-5', '6-7', '8-9')) +
    geom_smooth(col = 'grey50', method = 'lm') +   
    labs(x = 'Elevation range (km)', y = 'Relative HWI') +
    theme_classic() + 
  # add the stats
    ggtext::geom_richtext(aes(x=2.22, y=0.82, label = stat.c), label.color = NA, size = 5) +
    theme(legend.position = 'none', 
          axis.title = element_text(size = 18), axis.text = element_text(size = 16),
          plot.margin = margin(t=5, l=1, r=30))

q2= ggplot(datbin.plot, aes(elevbin, mid.wa)) + 
  # median HWA
    geom_point(aes(size = n), col = 'brown') +  
      # add dot size legend
      scale_size_continuous(breaks = c(1000, 2000, 3000), 
                            labels = c('= 1,000 species', '= 2,000 species', '= 3,000 species')) + 
  # 1st and 3rd quartiles of HWA
    geom_errorbar(aes(ymin = lower.wa, ymax = upper.wa), col = 'grey50') + 
      scale_y_continuous(limits = c(-0.32, 0.34), breaks = c(-0.2, 0, 0.2)) +
      scale_x_continuous(breaks = c(0,2,4,6,8), labels = c('0-1','2-3', '4-5', '6-7', '8-9')) +
    geom_smooth(col = 'grey50', method = 'lm') + 
    labs(x = 'Elevation range (km)', y = 'Relative HWA') +
    theme_classic() + 
  # add the stats
    ggtext::geom_richtext(aes(x=2, y=0.27, label = stat.d), label.color = NA, size = 5) +
    theme(legend.position = c(0.82, 0.22), legend.title = element_blank(),
          legend.text=element_text(size=12), legend.margin = margin(c(0,0,0,0)),
          axis.title = element_text(size = 18), axis.text = element_text(size = 16),
          plot.margin = margin(t=5, l=1, r=25))



## figure 2
pdf('Figures/Figure2.pdf', width = 11, height = 7.2)
top = plot_grid(t2, t4,  ncol = 2, scale=0.97,
                 labels =  c(' A', ' B'), label_size = 18)
                       
bottom = plot_grid(q1, q2, ncol = 2, scale = 0.95,
                   labels =  c(' C', ' D'), label_size = 18)

plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 0.75))
dev.off()



############## ......................... Figure 3 .........................##############

######## // Fig.3A-B: forest plots // ########

library(ggbreak)

# function for averaging global PGLS model results and plotting
fig3 <- function(mod_no1, mod_no2, col, title=NULL){
  
  # read the all species model results
  coef_full = readRDS(coefs[mod_no1]) 
    # calculate the average model results using Rubin's rule (Nakagawa & De Villemereuil 2019;
    # https://doi.org/10.1093/sysbio/syy089)
    pooled_full <- pool.table(coef_full, type = "all")
    # add model labels
    pooled_full$model =  paste0('1', coefs[mod_no1])
  
  # read the landbird model results
  coef_reduce = readRDS(coefs[mod_no2])
    pooled_reduce <- pool.table(coef_reduce, type = "all")
    pooled_reduce$model =  paste0('0.', coefs[mod_no2])
  
  pooled = rbind(pooled_full, pooled_reduce)
  
  ## add sample sizes
  N1 = scales::comma(pooled_full$dfcom[1] + length(pooled_full$term))
  N2 = scales::comma(pooled_reduce$dfcom[1] + length(pooled_reduce$term))
  
  ## plot
  pooled$sig = ifelse(pooled$conf.low *  pooled$conf.high > 0, 1, 0.3)
  pooled = filter(pooled, term != '(Intercept)')
  
  p = ggplot2::ggplot(pooled, aes(term, estimate, shape = model, alpha =  I(sig))) + 
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    linewidth = 1.2, col = col,
                    position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, lty = 'solid', col = 'grey85') + 
    theme_classic() + 
    scale_x_discrete(limits = rev(y_orders), labels = rev(y_labs)) +
    scale_shape_manual(values = c(17, 16),
                       labels = c(paste0('Non-migratory landbirds (n=', N2, ')'), 
                                  paste0('All species (n=', N1, ')')),
                       guide = guide_legend(reverse = TRUE,
                                            label.position = 'left',
                                            label.hjust = 1)) +
    labs(col = '', x = '', y = 'Standardised effect size') +
    ggtitle(title)
  
  return(p)
  
}


## re-order the variables
y_orders = c('Elev', 'Lat', 'TempVar',  'Habitat.Openness', 
             'log.Mass', 'Flight.Modesoaring', 'Aerial.Lifestyle', 'Migration',
             'Trophic.Level2ry consumer', 'Elev:Flight.Modesoaring')
## add labels
y_labs = c('Elevation', 'Latitude', 'Temp. seasonality', 'Habitat openness',
           'Body mass', 'Flight mode (soaring)', 'Aerial lifestyle', 'Migration', 
           'Trophic level (2ry consumer)', 'Elevation : Soaring')



setwd("X:/home/_Wings_/HPC/global_coefs100/")
coefs = list.files()


# HWI results (Fig.3A)
f3a = fig3(6,8, '#2837a8') +
  annotate(geom = "rect", xmin = 0.5, xmax = 6.45, ymin = -0.44, ymax = -0.173,
           fill = "#c4c1d4", colour = "white", alpha = 0.2) + 
  annotate(geom = "rect", xmin = 6.55, xmax = 10.5, ymin = -0.44, ymax = -0.173,
           fill = "lightblue", colour = "white", alpha = 0.23) +   
  coord_flip(ylim = c(-0.15, 0.25), clip = 'off') + 
  theme(legend.position = 'none', 
        axis.text = element_text(size = 14),  
        axis.title.x.bottom = element_text(size = 15, vjust = -1.2))
ggsave('../Figures/Fig. 3A.svg', f3a, width = 7.2, height = 4.1)


# HWA results (Fig.3B)
f3b = fig3(24,26, 'coral') + scale_y_continuous(limits = c(-0.1, 1.03),
                                 breaks = c(-0.1, 0, 0.1, 0.2)) + 
  # body size has large effect size compared to other predictors, break the axis for visualisation 
  scale_y_break(c(0.16, 0.86), ticklabels = c(0.9, 1), space = 0.3) +
  theme(axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.x.bottom = element_text(size = 14), 
        axis.title.x.bottom = element_text(size = 14, hjust = 0.2, vjust = 2), 
        legend.text = element_text(size = 10)) +
  coord_flip()
ggsave('../Figures/Fig. 3B.svg', f3b, width = 7.4, height = 4.1) 




######## table S1 (related to fig.3) ########

# for tidying up the model summary tables
num_tidy = function(table){
    # round the main estimates to 3 digits
    table[,c('estimate', 'std.error', 'conf.low', 'conf.high')] <- round(table[,c('estimate', 'std.error', 'conf.low', 'conf.high')], 3)
    # change p values to significance levels (or round to 2 digits if > 0.05) 
    for (i in 1:nrow(table)){
        if (!is.na(table$p.value[i])) {
               if (as.numeric(table$p.value[i]) < 0.001) {table$p.value[i] <- '< 0.001'} 
              try (if (between(as.numeric(table$p.value[i]), 0.001, 0.01)) {table$p.value[i] <- '< 0.01'} )
              try (if (between(as.numeric(table$p.value[i]), 0.01, 0.05)) {table$p.value[i] <- '< 0.05'} )
              try (if (as.numeric(table$p.value[i]) >= 0.05) {table$p.value[i] <- round(as.numeric(table$p.value[i]), 2)})
        }
    }
    colnames(table) <- c('Predictor', 'Estimate', 'Std. error', 'Lower 95% CI', 'Higher 95% CI', 'Pt value', 'Model')
    return(table)
}

# for producing the long supplementary tables 1 and 3
Supp_table <- function(a, b, c, d, e){
  
  ## all species (line 1)
    coef1 = readRDS(coefs[a]) %>% pool.table(type = "all")
    coef1$model =  paste0('5', coefs[a])  # add model labels
  ## non-migratory landbirds (line 2)
    coef2 = readRDS(coefs[b]) %>% pool.table(type = "all")
    coef2$model =  paste0('4.', coefs[b])
  ## high data certainty (line 3)
    coef3 = readRDS(coefs[c]) %>% pool.table(type = "all")
    coef3$model =  paste0('3.', coefs[c])
  ## mean elevation (line 4) 
    coef4 = readRDS(coefs[d]) %>% pool.table(type = "all")
    coef4$model =  paste0('2.', coefs[d])
  ## absolute wing changes (line 5) 
    coef5 = readRDS(coefs[e]) %>% pool.table(type = "all")
    # add body mass back to predictors (all values being NA) for table formatting
    coef5 = rbind(coef5, NA)
    levels(coef5$term) = c(levels(coef5$term), 'log.Mass')
    coef5$term[11] <- 'log.Mass'
    coef5$model =  paste0('1.', coefs[e])
   
    pooled = rbind(coef1, coef2, coef3, coef4, coef5)
    
    SI_table = pooled %>% 
      # re-order the variables
      mutate(term = factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term) %>% 
      select(c(term, estimate, std.error, conf.low, conf.high, p.value, model)) 
}


### 1. HWI results
Table_S1_HWI = Supp_table(6, 8, 3, 4, 1) %>% num_tidy() %>% write.csv('../SI tables/Table_S1_HWI.csv', row.names = F) # (ignore errors)

### 2. HWA results
Table_S1_HWA = Supp_table(24, 26, 21, 22, 19) %>% num_tidy() %>% write.csv('../SI tables/Table_S1_HWA.csv', row.names = F)

    ## (adding direct estimates for soaring species)
    HWI1 = readRDS(coefs[5]) %>% pool.table(type = "all") %>% filter(term == 'Elev')
    HWI2 = readRDS(coefs[7]) %>% pool.table(type = "all") %>% filter(term == 'Elev')  
    HWI5 = readRDS(coefs[2]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 
    
    HWA1 = readRDS(coefs[23]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 
    HWA2 = readRDS(coefs[25]) %>% pool.table(type = "all") %>% filter(term == 'Elev')  
    HWA5 = readRDS(coefs[20]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 

    

######## table S2 (related to fig.3) ########
    
### Body mass changes
    
## all species (line 1)
aMf = readRDS(coefs[9]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[9], term =  factor(term, levels = c('(Intercept)', y_orders[-5]), labels = c('(Intercept)', y_labs[-5]))) 

## non-migratory landbirds (line 2)
aMf_l = readRDS(coefs[10]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[10], term =  factor(term, levels = c('(Intercept)', y_orders[-5]), labels = c('(Intercept)', y_labs[-5]))) 

table2 = rbind(aMf, aMf_l) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod)) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI tables/Table_S2_mass.csv', row.names = F)



######## table S3 (related to fig.3) ########

### 1. Wing length results
Table_S3_WL = Supp_table(32, 34, 29, 30, 27) %>% num_tidy() %>% write.csv('../SI tables/Table_S3_WL.csv', row.names = F) # (ignore the errors)

### 2. Wing width results
Table_S3_SL = Supp_table(16, 18, 13, 14, 11) %>% num_tidy() %>% write.csv('../SI tables/Table_S3_SL.csv', row.names = F)

    ## (adding direct estimates for soaring species)
    WL1 = readRDS(coefs[31]) %>% pool.table(type = "all") %>% filter(term == 'Elev')
    WL2 = readRDS(coefs[33]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 
    WL5 = readRDS(coefs[28]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 
    
    SL1 = readRDS(coefs[15]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 
    SL2 = readRDS(coefs[17]) %>% pool.table(type = "all") %>% filter(term == 'Elev')  
    SL5 = readRDS(coefs[12]) %>% pool.table(type = "all") %>% filter(term == 'Elev') 



######## R2 and Pagel's lambda values ########
    
setwd('X:/home/_Wings_/HPC/raw_models_100_global/')

r2_lambda_table <- data.frame(mod=NA, run=NA, condi_R2=NA, lambda=NA)
for (i in 1:length(dir())){       
  mod_set_i = readRDS(dir()[i])
  for (j in 1:length(mod_set_i)){
    line_i = c(dir()[i], j, R2_lik(mod_set_i[[j]]), summary(mod_set_i[[j]])$optpar)
    r2_lambda_table = rbind(r2_lambda_table, line_i)
  }
  rm(mod_set_i)
}

# calculate the mean
summary_mods <- r2_lambda_table[-1,] %>% group_by(mod) %>% 
  summarise(mean_condi_R2 = round(mean(as.numeric(condi_R2)), 3),
            mean_lambda = round(mean(as.numeric(lambda)), 3))
write.csv(summary_mods, '../SI tables/r2.lambda summary.csv')




######## // Fig.3D-E: within-family comparisons // ########

# use the more recent Clements taxonomy for this analysis

{ # add Clements taxonomy info
  cw = read.csv('Data S1_crosswalk.csv')
  cw = unique(cw[,c('Species_BirdTree', 'Species_Clements', 'Family_Clements')]) %>% na.omit()
  dat = left_join(dat, cw, by='Species_BirdTree')
  dat = select(dat, -c('Species_Clements')) %>% unique() # remove duplicates
}

# remove irrelevant species 
dat = filter(dat, Migration != 3, Habitat.type != 'Marine')


seg_plot = function(wing_metric){
  
      # calculate wing data for group 1 [Low elevation + High latitude]
      d1 = subset(dat, dat$MAX <= 3000 & abs(dat$Min.latitude) >= 30 & dat$Migration %in% c(1,2))
      t1 = d1 %>% group_by(Family_Clements) %>% summarise(mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                          mid.wa = median(rela.WA, na.rm = TRUE),
                                                          n=n(), cate = '1')
      # calculate wing data for group 2 [High elevation + Low latitude]
      d2 = subset(dat, dat$MAX >= 4000 & abs(dat$Max.latitude) <= 30 & dat$Migration %in% c(1,2))
      t2 = d2 %>% group_by(Family_Clements) %>% summarise(mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                          mid.wa = median(rela.WA, na.rm = TRUE),
                                                          n=n(), cate = '2')
      # merge the 2 groups by family and compare
      dat_plot = merge(t1, t2, by = 'Family_Clements')
      # remove families with too small (< 3) sample sizes
      dat_plot = filter(dat_plot, n.x >= 3 & n.y >= 3)
      
  
  if (wing_metric == 'HWI'){
    # show stats
    print(t.test(dat_plot$mid.hwi.y, dat_plot$mid.hwi.x, paired = TRUE))
    
    # calculate if the wing metric increases or decreases from group 1 to 2.
    dat_plot$del = dat_plot$mid.hwi.y - dat_plot$mid.hwi.x
    # use solid line for increases, dashed lines for decreases
    dat_plot$direction = ifelse(dat_plot$del > 0, 'solid', 'longdash')
    
    p = ggplot(dat_plot, aes(group = Family_Clements)) + 
      geom_segment(aes(x = cate.x,   y = mid.hwi.x,  xend = cate.y, yend = mid.hwi.y,
                       linetype = I(direction))) +
      geom_point(aes(x = cate.x, y = mid.hwi.x), col = '#2837a8') + 
      geom_point(aes(x = cate.y, y = mid.hwi.y), col = '#2837a8') + 
      scale_x_discrete(breaks = c(1,2),labels = c('Low elevation\nHigh latitude', 'High elevation\nLow latitude')) +
      scale_y_continuous(limits = c(-0.8, 1.5)) +
      labs(x = '', y = 'Relative hand-wing index', col = 'Family', size = '# of species') + 
      theme_classic() + theme(legend.position = 'none', axis.text = element_text(size = 12), axis.title = element_text(size = 12))
  }
  
  if (wing_metric == 'HWA'){
    # show stats
    print(t.test(dat_plot$mid.wa.y, dat_plot$mid.wa.x, paired = TRUE))
    
    dat_plot$direction = ifelse(dat_plot$mid.hwi.y > dat_plot$mid.hwi.x, 'solid', 'longdash')
    
    p = ggplot(dat_plot, aes(group = Family_Clements)) + 
      geom_segment(aes(x = cate.x,   y = mid.wa.x, xend = cate.y, yend = mid.wa.y,
                       linetype = I(direction))) +
      geom_point(aes(x = cate.x, y = mid.wa.x), col = 'coral') +
      geom_point(aes(x = cate.y, y = mid.wa.y), col = 'coral') + 
      scale_x_discrete(breaks = c(1,2),labels = c('Low elevation\nHigh latitude', 'High elevation\nLow latitude')) +
      scale_y_continuous(limits = c(-0.4, 0.4)) +
      labs(x = '', y = 'Relative hand-wing area', col = 'Family', size = '# of species') + 
      theme_classic() +  theme(legend.position = 'none', axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
      annotate('text', x = 1.5,  y = 0.37, 
               label = paste0('Non-migratory landbirds\nn=', sum(dat_plot$n.x), ' species              n=', sum(dat_plot$n.y), ' species'))
    
  }
  return(p)
}


### figure 3D-E
pd = seg_plot(wing_metric = 'HWI') 
pe = seg_plot(wing_metric = 'HWA')
ggsave('Figures/Fig. 3DE.svg', pd + pe, width = 9.4, height = 3.25)



######## figure S3 (related to fig.3D-E) ########

## adapted from seg_plot function for showing additional info (families, sample sizes, stats)
seg_plot_SI = function(wing_metric){
  
    d1 = subset(dat, dat$MAX <= 3000 & abs(dat$Min.latitude) >= 30 & dat$Migration %in% c(1,2))
    t1 = d1 %>% group_by(Family_Clements) %>% summarise(mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                        mid.wa = median(rela.WA, na.rm = TRUE),
                                                        n=n(), cate = '1')
    
    d2 = subset(dat, dat$MAX >= 4000 & abs(dat$Max.latitude) <= 30 & dat$Migration %in% c(1,2))
    t2 = d2 %>% group_by(Family_Clements) %>% summarise(mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                        mid.wa = median(rela.WA, na.rm = TRUE),
                                                        n=n(), cate = '2')
    
    dat_plot = merge(t1, t2, by = 'Family_Clements')
    
    # remove families with too small (< 3) sample sizes
    dat_plot = filter(dat_plot, n.x >= 3 & n.y >= 3)
  
  if (wing_metric == 'HWI'){
    
    # for adding stats
    t <- t.test(dat_plot$mid.hwi.y, dat_plot$mid.hwi.x, paired = TRUE)
    
    dat_plot$direction = ifelse(dat_plot$mid.hwi.y > dat_plot$mid.hwi.x, 'solid', 'longdash')
    
    p = ggplot(dat_plot, aes(group = Family_Clements)) + 
      geom_segment(aes(x = cate.x,   y = mid.hwi.x, xend = cate.y, yend = mid.hwi.y,
                       linetype = I(direction))) +
      geom_point(aes(x = cate.x, y = mid.hwi.x, size = n.x, col = Family_Clements)) + 
      geom_point(aes(x = cate.y, y = mid.hwi.y, size = n.y, col = Family_Clements)) + 
      scale_x_discrete(breaks = c(1,2),labels = c('Low elevation\nHigh latitude', 'High elevation\nLow latitude')) +
      scale_y_continuous(limits = c(-0.9, 2)) +
      labs(x = '', y = 'Relative hand-wing index') + 
      annotate(geom = 'text', x = 1.5, y = 1.9, size = 4.5, 
               label = paste0('Paired t-test\n95% CI = [', round(t$conf.int[1], 3), ', ', round(t$conf.int[2], 3),
                              ']\nt = ', round(t$statistic, 2), ', P = ', round(t$p.value, 2))) +
      theme_classic() + theme(legend.position = 'none', axis.text = element_text(size = 14), axis.title = element_text(size = 16))
  }
  
  
  if (wing_metric == 'HWA'){
    
    # for adding stats
    t = t.test(dat_plot$mid.wa.y, dat_plot$mid.wa.x, paired = TRUE)
    
    dat_plot$direction = ifelse(dat_plot$mid.hwi.y > dat_plot$mid.hwi.x, 'solid', 'longdash')
    
    p = ggplot(dat_plot, aes(group = Family_Clements)) + 
      geom_segment(aes(x = cate.x,   y = mid.wa.x, xend = cate.y, yend = mid.wa.y,
                       linetype = I(direction))) +
      geom_point(aes(x = cate.x, y = mid.wa.x, size = n.x, col = Family_Clements)) + 
      geom_point(aes(x = cate.y, y = mid.wa.y, size = n.y, col = Family_Clements)) + 
      scale_x_discrete(breaks = c(1,2),labels = c('Low elevation\nHigh latitude', 'High elevation\nLow latitude')) +
      scale_y_continuous(limits = c(-0.4, 0.32)) +
      scale_size_continuous(breaks = c(10, 20, 30), labels = c('= 10 species', '= 20 species', '= 30 species')) +
      labs(x = '', y = 'Relative hand-wing area', col = 'Clade (family)', size = 'Sample size per family') + 
      annotate(geom = 'text', x = 1.5, y = 0.3, size = 4.5,
               label = paste0('Paired t-test\n95% CI = [', round(t$conf.int[1], 3), ', ', round(t$conf.int[2], 3),
                              ']\nt = ', round(t$statistic, 2), ', P < 0.01')) +
      theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
      guides(colour = guide_legend(override.aes = list(size=3), order = 1),
             size = guide_legend(order = 0))
  }
  return(p)
}


## figure S3
pi = seg_plot_SI(wing_metric = 'HWI') 
pa = seg_plot_SI(wing_metric = 'HWA') 
ps3 = cowplot::plot_grid(pi, pa, rel_widths = c(1, 1.54))
ggsave(file = 'Figures/fig.S3.svg', plot=ps3, width=9.6, height=6)




############## ......................... Figure 4 .........................##############

##### Step 1: average coefficients from the 100 trees using Rubin's rule 

ave_100_by_group <- function(coef_table){
  
    coef_mean_ntree=c()
    for (bin in unique(coef_table$elevbin) ){
        for (y in unique(coef_table$Y)) {
          
            # this contains the 100 model results per elevation band per wing metric
            coef_bin_y_i <- coef_table %>% filter(elevbin == bin, Y == y)
            
            # format table to be used by pool.table function
            coef_bin_y_i$term = coef_bin_y_i$predictor
            coef_bin_y_i$df.residual = coef_bin_y_i$sampleSize - length(unique(coef_bin_y_i$term))
            colnames(coef_bin_y_i)[c(1,2)] <- c('estimate', 'std.error')
            
            pooled_bin_y_i = pool.table(coef_bin_y_i, type = "all")
            
            # add back the model info
            pooled_bin_y_i$Y = y
            pooled_bin_y_i$elevbin = bin
            pooled_bin_y_i$sampleSize = unique(coef_bin_y_i$sampleSize)
            
            # combine results
            coef_mean_ntree = rbind(coef_mean_ntree, pooled_bin_y_i)
        }
    }
    return(coef_mean_ntree)
}



##### Step 2: plot the coefficients

plot_by_group <- function(coef){
  
  # use paler colours for non-significant results 
  coef$sig = ifelse(coef$p.value < 0.05, 1, 0.3)
  coef$elevbin = as.character(coef$elevbin)
  # select the main predictor to be shown
  coef_elev <<- filter(coef, term %in% c('Elev'))
  # re-order response variables
  coef_elev$Y <- factor(coef_elev$Y, levels = c('log.Sec1', 'log.Wing.Length', 'log.WA', 'HWI'))

  p = ggplot(coef_elev, aes(elevbin, estimate, col = Y, alpha = I(sig))) + 
          geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                          linewidth = 1.4, position = position_dodge(width = 0.55)) +
          geom_hline(yintercept = 0, lty = 'solid', col = 'grey80') + 
          labs(y = 'Effect size of elevation')+
          theme_classic() +
          scale_color_manual(name = '', guide = guide_legend(reverse = TRUE),
                             values=c('brown', 'olivedrab', 'coral', '#2837a8'),
                             labels = c('Wing width', 'Wing length','Hand-wing area','Hand-wing index')) +
          coord_flip() 
    
    return(p)
    
}


  
setwd("X:/home/_Wings_/HPC/elev_band_coef100/")

## figure 4
f4 = ave_100_by_group(readRDS('all species.rds')) %>% plot_by_group() +
    scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   # add sample size info for each elevation band
                   labels=c(paste0('0-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'),
                            paste0('5-8\n(n=', scales::comma(unique(coef_elev$sampleSize)[6]), ')'))) +
  theme(axis.text = element_text(size=14),  axis.title = element_text(size=16),
        legend.title = element_text(size=12), legend.text = element_text(size=12))
ggsave('../Figures/Fig. 4.svg', f4, width = 9.5, height = 5.3)



######## figure S4 (related to fig.4) ########

# (A) absolute wing changes
p1 = ave_100_by_group(readRDS('absolutue wing changes.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'),
                            paste0('5-8\n(n=', scales::comma(unique(coef_elev$sampleSize)[6]), ')'))) +
  theme(legend.position = 'none',  plot.margin = margin(l=7, t=7, r=0, b=30),
        axis.text = element_text(size=18), axis.title = element_text(size=24))


# (B) 4-km bandwith
p2 = ave_100_by_group(readRDS('bandwidth = 4km.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-8\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'))) +
  scale_y_continuous(breaks = c(0, 0.02, 0.04)) +
  theme(legend.position = 'none', plot.margin = margin(l=0, t=7, r=0, b=30),
        axis.text = element_text(size=18), axis.title = element_text(size=24))


# (C) 2-km bandwith
p3 = ave_100_by_group(readRDS('bandwidth = 2km.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-2\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'),
                            paste0('5-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[6]), ')'))) +
  theme(legend.position = c(0.7, 0.15), legend.title = element_blank(), legend.text = element_text(size=18.5),
        plot.margin = margin(l=0, t=7, r=7, b=30),
        axis.text = element_text(size=18),  axis.title = element_text(size=24))


# (D) landbird model
p4 = ave_100_by_group(readRDS('landbird model.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'),
                            paste0('5-8\n(n=', scales::comma(unique(coef_elev$sampleSize)[6]), ')'))) +
  theme(legend.position = 'none', plot.margin = margin(l=7, t=21, r=0),
        axis.text = element_text(size=18), axis.title = element_text(size=24))


# (E) high certainty species
p5 = ave_100_by_group(readRDS('high certainty.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'),
                            paste0('4-7\n(n=', scales::comma(unique(coef_elev$sampleSize)[5]), ')'),
                            paste0('5-8\n(n=', scales::comma(unique(coef_elev$sampleSize)[6]), ')'))) + 
  theme(legend.position = 'none', plot.margin = margin(l=0, t=21, r=0),
        axis.text = element_text(size=18), axis.title = element_text(size=24))


# (F) use mean elevation instead of max
p6 = ave_100_by_group(readRDS('mean.elevation.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', scales::comma(unique(coef_elev$sampleSize)[1]), ')'),
                            paste0('1-4\n(n=', scales::comma(unique(coef_elev$sampleSize)[2]), ')'),
                            paste0('2-5\n(n=', scales::comma(unique(coef_elev$sampleSize)[3]), ')'), 
                            paste0('3-6\n(n=', scales::comma(unique(coef_elev$sampleSize)[4]), ')'))) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1)) +
  theme(legend.position = c(0.7, 0.15), legend.title = element_blank(), legend.text = element_text(size=18.5),
        plot.margin = margin(l=0, t=21, r=7),
        axis.text = element_text(size=18), axis.title = element_text(size=24))


## figure s4
px6 = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, rel_widths = c(1.05, 0.97, 1), 
                        labels = c('A', 'B', 'C', 'D', 'E', 'F'), label_size = 28)
ggsave(file = '../Figures/fig.S4.svg', plot = px6, width = 16.2, height = 12.9)




######## figure S6 (related to fig.4) ########

dat = read.csv('../Data S1_BirdTree.csv') %>% rename(MAX = 'Max.Elevation.1', MEAN = 'Mean.Elevation.1', MIN = 'Min.Elevation.1')
dat = filter(dat, Flight.Mode != 'flightless', Migration != 3)

## assign elevation bands to each species (same as fig.4)
{
  dat$elevbin1 = ifelse(dat$MAX <= 3000, '0-3', 0)
  dat$elevbin2 = ifelse(dat$MAX >= 1000 & dat$MAX <= 4000, '1-4', 0)
  dat$elevbin3 = ifelse(dat$MAX >= 2000 & dat$MAX <= 5000, '2-5', 0)
  dat$elevbin4 = ifelse(dat$MAX >= 3000 & dat$MAX <= 6000, '3-6', 0)
  dat$elevbin5 = ifelse(dat$MAX >= 4000 & dat$MAX <= 7000, '4-7', 0)
  dat$elevbin6 = ifelse(dat$MAX >= 5000, '5-8', 0)

}

dat_long = rbind(dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin1')] %>% rename(elevbin = 'elevbin1'),
                 dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin2')] %>% rename(elevbin = 'elevbin2'),
                 dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin3')] %>% rename(elevbin = 'elevbin3'),
                 dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin4')] %>% rename(elevbin = 'elevbin4'),
                 dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin5')] %>% rename(elevbin = 'elevbin5'),
                 dat[,c('Species_BirdTree', 'MAX', 'MEAN', 'MIN', 'elevbin6')] %>% rename(elevbin = 'elevbin6'))
dat_long = filter(dat_long, elevbin != '0')

# most bands contain too many species to be plotted individually, randomly choose a subset of 500 species
# (500 is a lot so the subsets can still represent data distribution)
{
rows1 = which(dat_long$elevbin == '0-3') %>% sample(500)
rows2 = which(dat_long$elevbin == '1-4') %>% sample(500)
rows3 = which(dat_long$elevbin == '2-5') %>% sample(500)
rows4 = which(dat_long$elevbin == '3-6') %>% sample(500)
rows5 = which(dat_long$elevbin == '4-7') %>% sample(500)
rows6 = which(dat_long$elevbin == '5-8')

dat_long_sub = rbind(dat_long[rows1, ], dat_long[rows2, ], dat_long[rows3, ],
                     dat_long[rows4, ], dat_long[rows5, ], dat_long[rows6, ])

}

# use the full dataset for making boxplots
dat_longer_full = pivot_longer(dat_long, cols = 2:4, names_to = 'Variables', values_to = 'Elevation')
# use the subset for making scatter plots
dat_longer_sub = pivot_longer(dat_long_sub, cols = 2:4, names_to = 'Variables', values_to = 'Elevation')
  
# standardise units
dat_longer_full$Elevation = dat_longer_full$Elevation/1000
dat_longer_sub$Elevation = dat_longer_sub$Elevation/1000


Box_Mean <- function(x) {
  v <- c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

ggplot(data = dat_longer_full, aes(x=elevbin, y=Elevation, fill=Variables)) +
  # add boxplots using the full data
  stat_summary(fun.data = Box_Mean, geom = "boxplot", position = position_dodge(width = 1), width = 0.7, col = 'black') + 
  # add scatter plots using the subsets
  geom_point(data = dat_longer_sub, aes(col = Variables), alpha = 0.5, size = 0.5, col = 'grey70', 
             position = position_jitterdodge(jitter.width = 0.18, jitter.height = 0.02, dodge.width = 1) ) +
  labs(fill = '', x = 'Elevation band (km)', y = 'Elevation (km)') +
  scale_fill_manual(values = c('#0096c7', '#8ecae6','#caf0f8' ), 
                    labels = c('Maximum elevation', 'Mean elevation', 'Minimum elevation'))+
  scale_y_continuous(limits = c(-0.500, 8.500)) +
  # add sample sizes
  annotate(geom = 'text', x = 3.5, y = -0.450, size = 5.1,
           label = paste('n =', scales::comma(table(dat_long$elevbin)[1]), '        ',
                         'n =', scales::comma(table(dat_long$elevbin)[2]), '        ',
                         'n =', scales::comma(table(dat_long$elevbin)[3]), '         ',
                         'n =', scales::comma(table(dat_long$elevbin)[4]), '         ',
                         'n =', scales::comma(table(dat_long$elevbin)[5]), '         ',
                         'n =', scales::comma(table(dat_long$elevbin)[6]))) +
  theme_classic() +
  theme(axis.title = element_text(size = 17), axis.text = element_text(size = 15),
        legend.position = c(0.18, 0.85), legend.text = element_text(size = 13),
        legend.key.size = unit(1, 'cm'))

ggsave('../Figures/fig.S6.svg', width = 9.5, height = 7)




############## ................... Other supplementary figures ...................##############

setwd('X:/home/_Wings_/HPC/')

########  Fig. S1 ######## 

dat = read.csv('Data S1_BirdTree.csv') %>% rename(MAX = 'Max.Elevation.1', MEAN = 'Mean.Elevation.1')
{# add order-level taxonomic info
  cw = read.csv('Data S1_crosswalk.csv')
  cw = unique(cw[,c('Species_BirdTree', 'Order_Clements')]) %>% na.omit()
  dat = left_join(dat, cw, by='Species_BirdTree')
}

# group species into 10-degree latitudinal bins
dat$lat_bin = abs(dat$Centroid.latitude) %/% 10
# select high-elevation species only
dat=filter(dat, MAX >= 4000)

## figure S1A
non_passer = dat %>% filter(Order_Clements != 'Passeriformes') %>% select(all_of(c('Migration', 'lat_bin')))
P1 = ggplot(non_passer, aes(x=lat_bin, fill=as.character(Migration))) + 
  scale_x_continuous(breaks = seq(0, 16) - 0.5, labels = seq(0, 16) * 10) +
  # show the % of migrants per latitudinal bin
  geom_bar(position = 'fill') + scale_fill_manual(values = c('#deebf7', '#9ec9e0', '#3d77a1')) +
  labs(fill = 'Migration', x = 'Latitude', y = 'Proportion') +
  theme_classic() + theme(legend.position = 'none', axis.title = element_text(size = 12)) +
  # add sample sizes
  annotate(geom = 'text', x = 3.4, y = 1.05, col = 'grey50', size = 3.3, 
           label = paste(table(non_passer$lat_bin)[1],'   ',
                           table(non_passer$lat_bin)[2],'   ', 
                           table(non_passer$lat_bin)[3],'   ', 
                           table(non_passer$lat_bin)[4],'   ',
                           table(non_passer$lat_bin)[5],'   ',
                           table(non_passer$lat_bin)[6],'    ', 
                           table(non_passer$lat_bin)[7],'     ',
                           table(non_passer$lat_bin)[8]))
## figure S1B
passer = dat %>% filter(Order_Clements == 'Passeriformes') %>% select(all_of(c('Migration', 'lat_bin')))
# no passerines has a centroid latitude > 70 degrees, add an extra column here to directly compare with non-passerines
passer = rbind(passer, c(NA, 7))  
P2 = ggplot(passer, aes(x=lat_bin, fill=as.character(Migration))) + 
  scale_x_continuous(breaks = seq(0, 16) - 0.5, labels = seq(0, 16) * 10) +
  # show the % of migrants per latitudinal bin
  geom_bar(position = 'fill') + scale_fill_manual(values = c('#deebf7', '#9ec9e0', '#3d77a1')) +
  labs(fill = 'Migration', x = 'Latitude', y = '') +
  theme_classic() + theme(axis.title = element_text(size = 12)) +
  # add sample sizes
  annotate(geom = 'text', x = 3.45, y = 1.05, col = 'grey50', size = 3.3, 
           label = paste(table(passer$lat_bin)[1],'   ',
                         table(passer$lat_bin)[2],'   ', 
                         table(passer$lat_bin)[3],'  ', 
                         table(passer$lat_bin)[4],'  ',
                         table(passer$lat_bin)[5],'   ',
                         table(passer$lat_bin)[6],'   ', 
                         table(passer$lat_bin)[7],'    ',
                         '0'))

p = cowplot::plot_grid(P1, P2, rel_widths = c(1, 1.25))
ggsave(file="Figures/fig.S1.svg", plot=p, width=7.7, height=3.2)




########  Fig. S2 ######## 

dat = read.csv('Data S1_BirdTree.csv') %>% rename(MAX = 'Max.Elevation.1', MEAN = 'Mean.Elevation.1')

### fig.S2C (HWI validation)

dat_shape = na.omit(dat[,c('Hand.Wing.Index','Wing.aspect.ratio')])

cor1 = cor.test(dat_shape$Hand.Wing.Index, dat_shape$Wing.aspect.ratio)

ggplot(dat_shape, aes(Wing.aspect.ratio, Hand.Wing.Index)) + theme_classic()+
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18)) +
  labs(x='Aspect ratio', y = 'Hand-wing index') +
  annotate('text', label=paste0("Pearson's correlation = ", round(cor1$estimate, 2),
                                '\np < 0.001 \nn = ', nrow(dat_shape), ' species'), 
           hjust=0, x = 8.4, y = 12, size = 5.5)
ggsave('Figures/fig.S2C.svg', width = 5, height = 5)


### fig.S2D (HWA validation)

dat_area = na.omit(dat[,c('Total.wing.area','Hand.Wing.Area')]) %>% log()

cor2 = cor.test(dat_area$Total.wing.area, dat_area$Hand.Wing.Area)

ggplot(dat_area, aes(Total.wing.area, Hand.Wing.Area)) + theme_classic()+
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18)) +
  labs(x='Total wing area (log)', y = 'Hand-wing area (log)') +
  annotate('text', label=paste0("Pearson's correlation = ", round(cor2$estimate, 2),
                                '\np < 0.001 \nn = ', nrow(dat_area), ' species'), 
           hjust=0, x = 9.5, y = 7.5, size = 5.5)

ggsave('Figures/fig.S2D.svg', width = 5, height = 5)




########  Fig. S5 ######## 

elev = read.csv('Raw elevation data 3 sources.csv')

##### top panels: comparing raw elevation data

# get elevation data that consider all ranges of a species
elev$all_ranges_max = ifelse(is.na(elev$new_max), elev$WB_max, elev$new_max)
elev$all_ranges_min = ifelse(is.na(elev$new_min), elev$WB_min, elev$new_min)
elev$all_ranges_mean = (elev$all_ranges_max + elev$all_ranges_min)/2
# get elevation data that consider breeding ranges only
elev$QJ_mean.breeding. = (elev$QJ_max.breeding. + elev$QJ_min.breeding.)/2

# fig. S5A (max elevation)
elev_plot = na.omit(elev[, c('all_ranges_max', 'QJ_max.breeding.', 'Migration')])
elev_plot$Migration = ifelse(elev_plot$Migration == 3, 'Yes', 'No')
cor2 = cor.test(elev_plot$all_ranges_max, elev_plot$QJ_max.breeding.)
pa=ggplot(elev_plot, aes(all_ranges_max/1000, QJ_max.breeding./1000, col=Migration)) + theme_classic()+
      geom_point(size = 1, alpha = 0.5) +
      theme(axis.title = element_text(size=20),
            axis.text = element_text(size=18), legend.position = 'none') +
      labs(x='Elevation - all ranges (km)', y = 'Elevation - breeding range (km)') +
      scale_x_continuous(limits = c(0, 8.5)) +
      scale_color_manual(values = c('Black', 'red')) +
      annotate('text', label=paste0("Pearson's correlation = ", round(cor2$estimate, 2),
                                    '\np < 0.001 \nn = ', scales::comma(nrow(elev_plot)), ' species'), 
               hjust=0, x = 0, y = 7.7, size = 6)

# fig. S5B (mean elevation)
elev_plot = na.omit(elev[, c('all_ranges_mean', 'QJ_mean.breeding.', 'Migration')])
elev_plot$Migration = ifelse(elev_plot$Migration == 3, 'Yes', 'No')
cor3 = cor.test(elev_plot$all_ranges_mean, elev_plot$QJ_mean.breeding.)
pb=ggplot(elev_plot, aes(all_ranges_mean/1000, QJ_mean.breeding./1000, col=Migration)) + theme_classic()+
      geom_point(size = 1, alpha = 0.5) +
      theme(axis.title = element_text(size=20),
            axis.text = element_text(size=18), legend.position = 'none') +
      labs(x='Elevation - all ranges (km)', y = 'Elevation - breeding range (km)') +
      scale_x_continuous(limits = c(0, 6.2)) +
      scale_y_continuous(limits = c(0, 6.2)) +
      scale_color_manual(values = c('Black', 'red')) +
      annotate('text', label=paste0("Pearson's correlation = ", round(cor3$estimate, 2),
                                    '\np < 0.001 \nn = ', scales::comma(nrow(elev_plot)), ' species'), 
               hjust=0, x = 0, y = 5.8, size = 6)


##### bottom panels: comparing merged elevation data

dat = read.csv('Data S1_BirdTree.csv') 

# fig. S5C (max elevation)
elev_plot = na.omit(dat[, c('Max.Elevation.1', 'Max.Elevation.2')])
cor4 = cor.test(elev_plot$Max.Elevation.1, elev_plot$Max.Elevation.2)
pc=ggplot(elev_plot, aes(Max.Elevation.1/1000, Max.Elevation.2/1000)) + theme_classic()+
      geom_point(size = 1, alpha = 0.5) +
      theme(axis.title = element_text(size=20), axis.text = element_text(size=18), plot.margin = margin(t=25)) +
      labs(x='Maximum elevation 1 (km)', y = 'Maximum elevation 2 (km)') +
      scale_x_continuous(limits = c(0, 8.5)) +
      annotate('text', label=paste0("Pearson's correlation = ", round(cor4$estimate, 2),
                                    '\np < 0.001 \nn = ', scales::comma(nrow(elev_plot)), ' species'), 
               hjust=0, x = 0, y = 7.7, size = 6)

# fig. S5D (mean elevation)
elev_plot = na.omit(dat[, c('Mean.Elevation.1', 'Mean.Elevation.2')])
cor5 = cor.test(elev_plot$Mean.Elevation.1, elev_plot$Mean.Elevation.2)
pd=ggplot(elev_plot, aes(Mean.Elevation.1/1000, Mean.Elevation.2/1000)) + theme_classic()+
      geom_point(size = 1, alpha = 0.5) +
      theme(axis.title = element_text(size=20), axis.text = element_text(size=18), plot.margin = margin(t=25)) +
      labs(x='Mean elevation 1 (km)', y = 'Mean elevation 2 (km)') +
      scale_x_continuous(limits = c(0, 6.2)) +
      scale_y_continuous(limits = c(0, 6.2)) +
      annotate('text', label=paste0("Pearson's correlation = ", round(cor5$estimate, 2),
                                    '\np < 0.001 \nn = ', scales::comma(nrow(elev_plot)), ' species'), 
               hjust=0,  x = 0, y = 5.8, size = 6)


## figure S5
ps5 = cowplot::plot_grid(pa, pb, pc, pd, ncol = 2, scale = 0.97, labels = c('', '', 'C', 'D'), label_size = 25)
ggsave('Figures/fig.S5.svg', plot = ps5, width = 11, height = 10.5)

