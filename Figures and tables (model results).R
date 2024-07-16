setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/data")

library(tidyverse)
library(ggtext)
library(cowplot)


{
library(ape)
  library(caper)
library(phylolm)
  library(mice)
  library(rr2)
}


########## ................. Figure. 2 Raw patterns #########################

dat = read.csv('Analy data main.csv')

# calculate relative wing metrics for plotting
dat$rela.HWI = lm(dat$HWI ~ dat$log.Mass)$residuals
dat$rela.WA = lm(dat$log.WA ~ dat$log.Mass)$residuals


## panel a-b: scatter + density plot ####

cus_cols = c(# colour for the main scatter plot
             '#bedee8', '#012a7a', 
             # use slightly darker shades for the regression lines so they're more visible
             '#669eb0', '#153a85',
             # matching shades for Stn. Error bands
             '#669eb0', '#435d99',
             # density plots overlapping so must use transparent colours, using darker shades here.
             '#53d3fc', '#020a61')


p1=ggplot(dat, aes(MAX, rela.HWI, col=Flight_mode)) + geom_point(alpha=0.4) +
    theme_classic() +
    scale_colour_manual(values = c(cus_cols[1], cus_cols[2])) + 
  # set the x axis limits and expand here to align with side panels later.
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000), labels = c(0, 2, 4, 6, 8), limits = c(-400, 8700), expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
    labs(x = 'Maximum elevation (km)', y = 'Relative hand-wing index', col = 'Flight mode') +
  # add three regression lines:
    geom_smooth(data = dat[dat$Flight_mode == 'soar', ], col=cus_cols[4], fill=cus_cols[6], method = 'gam') +
    # remove marine soaring species  
    geom_smooth(data = dat[dat$Flight_mode == 'soar' & dat$Habitat != 'Marine', ], col=cus_cols[4], se=F, method = 'gam', linetype = 2, linewidth = 0.7) +
    geom_smooth(data = dat[dat$Flight_mode == 'flap', ], col=cus_cols[3], fill=cus_cols[5], method = 'gam') +
      theme(legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16))
           

p2=ggplot(dat, aes(MAX, rela.WA, col=Flight_mode)) + geom_point(alpha=0.4) +
    theme_classic() +
    scale_colour_manual(values = c(cus_cols[1], cus_cols[2])) + 
  # set the x axis limits and expand here to align with side panels later.
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000), labels = c(0, 2, 4, 6, 8), limits = c(-400, 8700), expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(x = 'Maximum elevation (km)', y = 'Relative hand-wing area', col = 'Flight mode') +
  # add three regression lines:
    geom_smooth(data = dat[dat$Flight_mode == 'soar', ], col=cus_cols[4], fill=cus_cols[6], method = 'gam') +
    # remove marine soaring species  
    geom_smooth(data = dat[dat$Flight_mode == 'soar' & dat$Habitat != 'Marine', ], col=cus_cols[4], se=F, method = 'gam', linetype = 2, linewidth = 0.7) +
    geom_smooth(data = dat[dat$Flight_mode == 'flap', ], col=cus_cols[3], fill=cus_cols[5], method = 'gam') +
      theme(legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16))


# make side density plots for panel a.

xdens1 <- axis_canvas(p1, axis = "x") +
    geom_density(data = dat, aes(x = MAX, fill = Flight_mode),
                alpha = 0.3, size = 0.25) +
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    scale_x_continuous(limits = c(-400, 8700))   # to align the x-axis of side panels to the main plot 

ydens1 <- axis_canvas(p1, axis = "y", coord_flip = TRUE) +
    geom_density(data = dat, aes(x = rela.HWI, fill = Flight_mode),
                  alpha = 0.3, size = 0.25)+
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    coord_flip()

t1 <- insert_xaxis_grob(p1, xdens1, grid::unit(.1, "null"), position = "top")
t2<- insert_yaxis_grob(t1, ydens1, grid::unit(.1, "null"), position = "right")


# make side density plots for panel b.
xdens2 <- axis_canvas(p2, axis = "x") +
    geom_density(data = dat, aes(x = MAX, fill = Flight_mode),
                alpha = 0.3, size = 0.25) +
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    scale_x_continuous(limits = c(-400, 8700))  # to align the x-axis of side panels to the main plot

ydens2 <- axis_canvas(p2, axis = "y", coord_flip = TRUE) +
    geom_density(data = dat, aes(x = rela.WA, fill = Flight_mode),
                  alpha = 0.3, size = 0.25)+
    scale_fill_manual(values = c(cus_cols[7], cus_cols[8])) +
    coord_flip()

t3 <- insert_xaxis_grob(p2, xdens2, grid::unit(.1, "null"), position = "top")
t4<- insert_yaxis_grob(t3, ydens2, grid::unit(.1, "null"), position = "right")


# get the legend
legends = get_legend(axis_canvas(p1, axis = "x")+
          geom_density(data = dat, aes(x = MAX, fill = Flight_mode),
                      alpha = 0.9, size = 0.2) +
          scale_fill_manual(values = c('#bedee8', '#153a85'),
                            labels = c('Flapping', 'Soaring')) +
          theme(legend.position = 'right') + theme_classic())




##  panel c-d: quantile lineplot  #####

datbin = dat
# split species into 1-km wide elevation bands based on their max elevation
datbin$elevbin = datbin$MAX %/% 1000  
# for each band, calculate the 1st, 2nd and 3rd quartiles of HWI and HWA
datbin.plot_full <- datbin %>% group_by(elevbin) %>% summarise( mid.hwi = median(rela.HWI, na.rm = TRUE),
                                                          lower.hwi = quantile(rela.HWI, 0.25, na.rm=TRUE),
                                                          upper.hwi = quantile(rela.HWI, 0.75, na.rm=TRUE),
                                                          mid.wa = median(rela.WA, na.rm = TRUE),
                                                          lower.wa = quantile(rela.WA, 0.25, na.rm=TRUE),
                                                          upper.wa = quantile(rela.WA, 0.75, na.rm=TRUE),
                                                          n=n())

datbin.plot = datbin.plot_full

# get the regression coefficients
stat.c = paste0('Slope = ',  round(lm(data = datbin.plot, mid.hwi ~ elevbin)$coefficients[2], 3), ', ', '*p*', ' < 0.001',
                          '<br> *R*<sup>2</sup> = ',  round(summary(lm(data = datbin.plot, mid.hwi ~ elevbin))$adj.r.squared, 2))

stat.d = paste0('Slope = ',  round(lm(data = datbin.plot, mid.wa ~ elevbin)$coefficients[2], 3), ', ', '*p*', ' < 0.01',
                          '<br> *R*<sup>2</sup> = ',  round(summary(lm(data = datbin.plot, mid.wa ~ elevbin))$adj.r.squared, 2))
  

q1= ggplot(datbin.plot, aes(elevbin, mid.hwi)) + 
    geom_point(aes(size = n), col = 'brown') +
    geom_errorbar(aes(ymin = lower.hwi, ymax = upper.hwi), col = 'grey50') + 
      scale_y_continuous(limits = c(-0.7, 1), breaks = c(-0.5, 0, 0.5)) +
      scale_x_continuous(breaks = c(0,2,4,6,8), labels = c('0-1','2-3', '4-5', '6-7', '8-9')) +
    geom_smooth(col = 'grey50', method = 'lm') +   
    labs(x = 'Elevation range (km)', y = 'Relative HWI') +
  theme_classic() + 
  ggtext::geom_richtext(aes(x=2.22, y=0.82, label = stat.c), label.color = NA, size = 5) +
  theme(legend.position = 'none', 
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        plot.margin = margin(t=5, r=30))

q2= ggplot(datbin.plot, aes(elevbin, mid.wa)) + 
    geom_point(aes(size = n), col = 'brown') +  
      scale_size_continuous(breaks = c(1000, 2000, 3000), 
                            labels = c('= 1000 species', '= 2000 species', '= 3000 species')) + 
    geom_errorbar(aes(ymin = lower.wa, ymax = upper.wa), col = 'grey50') + 
      scale_y_continuous(limits = c(-0.32, 0.34), breaks = c(-0.2, 0, 0.2)) +
      scale_x_continuous(breaks = c(0,2,4,6,8), labels = c('0-1','2-3', '4-5', '6-7', '8-9')) +
    geom_smooth(col = 'grey50', method = 'lm') + 
    labs(x = 'Elevation range (km)', y = 'Relative HWA') +
    theme_classic() + 
  ggtext::geom_richtext(aes(x=2, y=0.27, label = stat.d), label.color = NA, size = 5) +
  theme(legend.position = c(0.82, 0.22), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.margin = margin(c(0,0,0,0)),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        plot.margin = margin(t=5, r=25))



## combine all panels
pdf('../figures/Figure 2.pdf', width = 11, height = 7.2)
top = plot_grid(t2, t4,  ncol = 2, scale=0.97,
                 labels =  c(' a', ' b'), label_size = 18)
                       
bottom = plot_grid(q1, q2, ncol = 2, scale = 0.95,
                   labels =  c(' c', ' d'), label_size = 18)

plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 0.75))
dev.off()





## ..............  Figure 3. Global PGLS forest plots  ##################################

library(mice)
library(tidyverse)
library(ggbreak)

fig3 <- function(mod_no1, mod_no2, col, title=NULL){
  
  # calculate the average model results using Rubin's rules (Nakagawa & De Villemereuil 2019;
  # https://doi.org/10.1093/sysbio/syy089)
  coef_full = readRDS(coefs[mod_no1]) 
  pooled_full <- pool.table(coef_full, type = "all")
  pooled_full$model =  paste0('1', coefs[mod_no1])
  
  coef_reduce = readRDS(coefs[mod_no2])
  pooled_reduce <- pool.table(coef_reduce, type = "all")
  pooled_reduce$model =  paste0('0.', coefs[mod_no2])
  
  pooled = rbind(pooled_full, pooled_reduce)   
  ## show the results
  ## print(pooled)
  
  pooled$sig = ifelse(pooled$conf.low *  pooled$conf.high > 0, 1, 0.3)
  pooled = filter(pooled, term != '(Intercept)')
  
  p = ggplot(pooled, aes(term, estimate, shape = model, alpha =  I(sig))) + 
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    linewidth = 1.2, col = col,
                    position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, lty = 'solid', col = 'grey85') + 
    theme_classic() + 
    scale_x_discrete(limits = rev(y_orders),
                     labels = rev(y_labs)) +
    scale_shape_manual(values = c(17, 16),
                       labels = c(paste0('Non-migratory landbirds (n=', pooled_reduce$dfcom[1] + length(pooled_reduce$term), ')'), 
                                  paste0('All species (n=', pooled_full$dfcom[1] + length(pooled_full$term), ')')),
                       guide = guide_legend(reverse = TRUE,
                                            label.position = 'left',
                                            label.hjust = 1)) +
    labs(col = '', x = '', y = 'Standardised effect size') +
    ggtitle(title)
  
  return(p)
  
}


# re-order predictors and add labels
y_orders = c('Elev', 'Lat', 'TempVar',  'Habitat.Openness1.open', 
             'log.Mass', 'Flight_modesoar', 'AL.index1.frequent', 'Migration1.migratory',
             'Trophic.Level2ry consumer', 'Elev:Flight_modesoar')
y_labs = c('Elevation', 'Latitude', 'Temp. seasonality', 'Habitat openness',
           'Body mass', 'Flight mode (soaring)', 'Aerial lifestyle', 'Migration', 
           'Trophic level (2ry consumer)', 'Elevation : Soaring')



setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/global_coefs100/")
coefs = list.files()

# HWI results
pdf('C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/figures/figure 3a.pdf', width = 5.5, height = 4.45)
fig3(1, 2, '#192370') +
  annotate(geom = "rect", xmin = 0.5, xmax = 6.45, ymin = -0.485, ymax = -0.127,
           fill = "#c4c1d4", colour = "white", alpha = 0.2) + 
  annotate(geom = "rect", xmin = 6.55, xmax = 10.5, ymin = -0.485, ymax = -0.127,
           fill = "lightblue", colour = "white", alpha = 0.23) +   ## #99b84f, #8282b3, #73b0ff
  coord_flip(ylim = c(-0.1, 0.35), 
             clip = 'off') + 
  theme(legend.position = 'none', 
        axis.text = element_text(size = 12),  
        axis.title.x.bottom = element_text(size = 13.5, vjust = -1.2))
dev.off()


# HWA results
pdf('C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/figures/figure 3b.pdf', width = 6.8, height = 4.45)
fig3(9, 10, 'coral') + scale_y_continuous(limits = c(-0.1, 1.05),
                                 breaks = c(-0.1, 0, 0.1, 0.2)) + 
  scale_y_break(c(0.22, 0.85), ticklabels = c(0.9, 1), space = 0.3) +
  theme(axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.x.bottom = element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 14, hjust = 0.2, vjust = 2),
        legend.text = element_text(size = 10)) +
  coord_flip()
dev.off()



################# Supplementary analyses to Figure 3 ########################


### Figure S5. wing length (WL) and wing width (SL) results ##########

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/global_coefs100/")
coefs = list.files()

## WL
p1 = fig3(13, 14, 'olivedrab') + scale_y_continuous(limits = c(-0.1, 1.05),
                                 breaks = c(-0.1, 0, 0.1, 0.2)) + 
  scale_y_break(c(0.22, 0.85), ticklabels = c(0.9, 1), space = 0.3) +
  theme(axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.x.bottom = element_text(size = 12), 
        axis.text.y.left = element_text(size = 15), 
        axis.title.x.bottom = element_text(size = 15, hjust = 0.8, vjust = 2),
        legend.position = 'none') +
  coord_flip()


p2 = fig3(5, 6, 'brown') + scale_y_continuous(limits = c(-0.1, 1.05),
                                 breaks = c(-0.1, 0, 0.1, 0.2)) + 
  scale_y_break(c(0.22, 0.85), ticklabels = c(0.9, 1), space = 0.3) +
  theme(axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.x.bottom = element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 15, hjust = 0.25, vjust = 2)) +
  coord_flip()


pdf('C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/figures/SI/Fig. S5 WL+SL.pdf', width = 14, height = 5)
p1 + p2
dev.off()



###....................  Supplementary tables ............................###


num_tidy = function(table){
  table[,c('estimate', 'std.error', 'conf.low', 'conf.high')] <- round(table[,c('estimate', 'std.error', 'conf.low', 'conf.high')], 3)
  for (i in 1:nrow(table)){
    if (as.numeric(table$p.value[i]) < 0.001) {table$p.value[i] <- '< 0.001'} 
    try (if (between(as.numeric(table$p.value[i]), 0.001, 0.01)) {table$p.value[i] <- '< 0.01'} )
    try (if (between(as.numeric(table$p.value[i]), 0.01, 0.05)) {table$p.value[i] <- '< 0.05'} )
    try (if (as.numeric(table$p.value[i]) >= 0.05) {table$p.value[i] <- round(as.numeric(table$p.value[i]), 2)})
  }
  colnames(table) <- c('Predictor', 'Estimate', 'Std. error', 'Lower 95% CI', 'Higher 95% CI', 'Pt value', 'Model')
  return(table)
}


################ Table S1-2. data underlying figure 3 ####################

## Table S1. HWI results (all-species model, landbird model + one line for elevation (soar) estimate)

HWIf = readRDS(coefs[1]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[1], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
HWIs = readRDS(coefs[3]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[3]) # ref = soar
table_1a = rbind(HWIf, HWIs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

HWIf_l = readRDS(coefs[2]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[2], term =  factor(term, levels = c('(Intercept)', y_orders[-8]), labels = c('(Intercept)', y_labs[-8]))) %>% arrange(term)
HWIs_l = readRDS(coefs[4]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[4]) # ref = soar
table_1b = rbind(HWIf_l, HWIs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table1 = rbind(table_1a, table_1b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S1 HWI results (Related to figure 3).csv', row.names = F)  ## ignore the errors


## Table S2. HWA results (all-species model, landbird model + one line for elevation (soar) estimate)

WAf = readRDS(coefs[9]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[9], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
WAs = readRDS(coefs[11]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[11]) # ref = soar
table_2a = rbind(WAf, WAs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

WAf_l = readRDS(coefs[10]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[10], term =  factor(term, levels = c('(Intercept)', y_orders[-8]), labels = c('(Intercept)', y_labs[-8]))) %>% arrange(term)
WAs_l = readRDS(coefs[12]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[12]) # ref = soar
table_2b = rbind(WAf_l, WAs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table2 = rbind(table_2a, table_2b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S2 wing area results (Related to figure 3).csv', row.names = F)




####### Table S3-6: absolute wing changes (without size correction) #####################

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/global_coefs100__abs/")
coefs = list.files()

## Table S3. abs HWI results (all-species model, landbird model + one line for elevation (soar) estimate)

HWIf = readRDS(coefs[1]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[1], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
HWIs = readRDS(coefs[3]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[3]) # ref = soar
table_1a = rbind(HWIf, HWIs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

HWIf_l = readRDS(coefs[2]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[2], term =  factor(term, levels = c('(Intercept)', y_orders[-7]), labels = c('(Intercept)', y_labs[-7]))) %>% arrange(term)
HWIs_l = readRDS(coefs[4]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[4]) # ref = soar
table_1b = rbind(HWIf_l, HWIs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table1 = rbind(table_1a, table_1b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S3 absolute HWI changes (Related to figure 3).csv', row.names = F)  ## ignore the errors


## Table S4. abs HWA results (all-species model, landbird model + one line for elevation (soar) estimate)

WAf = readRDS(coefs[9]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[9], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
WAs = readRDS(coefs[11]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[11]) # ref = soar
table_2a = rbind(WAf, WAs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

WAf_l = readRDS(coefs[10]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[10], term =  factor(term, levels = c('(Intercept)', y_orders[-7]), labels = c('(Intercept)', y_labs[-7]))) %>% arrange(term)
WAs_l = readRDS(coefs[12]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[12]) # ref = soar
table_2b = rbind(WAf_l, WAs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table2 = rbind(table_2a, table_2b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S4 absolute wing area changes (Related to figure 3).csv', row.names = F)


## Table S5. abs wing length results (all-species model, landbird model + one line for elevation (soar) estimate)

HWIf = readRDS(coefs[13]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[13], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
HWIs = readRDS(coefs[15]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[15]) # ref = soar
table_1a = rbind(HWIf, HWIs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

HWIf_l = readRDS(coefs[14]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[14], term =  factor(term, levels = c('(Intercept)', y_orders[-7]), labels = c('(Intercept)', y_labs[-7]))) %>% arrange(term)
HWIs_l = readRDS(coefs[16]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[16]) # ref = soar
table_1b = rbind(HWIf_l, HWIs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table1 = rbind(table_1a, table_1b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S5 absolute wing length changes (Related to figure 3).csv', row.names = F)  ## ignore the errors


## Table S6. abs Wing width results (all-species model, landbird model + one line for elevation (soar) estimate)

WAf = readRDS(coefs[5]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[5], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
WAs = readRDS(coefs[7]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[7]) # ref = soar
table_2a = rbind(WAf, WAs[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

WAf_l = readRDS(coefs[6]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[6], term =  factor(term, levels = c('(Intercept)', y_orders[-7]), labels = c('(Intercept)', y_labs[-7]))) %>% arrange(term)
WAs_l = readRDS(coefs[8]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[8]) # ref = soar
table_2b = rbind(WAf_l, WAs_l[2,]) %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod))

table2 = rbind(table_2a, table_2b) %>% arrange(term) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S6 absolute wing width changes (Related to figure 3).csv', row.names = F)



############# Table S7-10. additional sensitivity analyses ###################

## Table S7-8. Additional HWI results (using high-certainty data; using mean elevation)

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/global_coefs100_SI/")
coefs = list.files()


HWI_high_cert = readRDS(coefs[1]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[1], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
table3 = HWI_high_cert %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod)) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S7 HWI results (Related to figure 3; supp models - high AL cert).csv', row.names = F)  ## ignore the errors


HWIf_mean_elev = readRDS(coefs[2]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[2], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)

table4 = HWIf_mean_elev %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod)) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S8 HWI results (Related to figure 3; supp models - mean elev).csv', row.names = F)  ## ignore the errors


## Table S9-10. Additional HWA results (using high-certainty data; using mean elevation)

WA_high_cert = readRDS(coefs[9]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[9], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)
table5 = WA_high_cert %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod)) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S9 HWA results (Related to figure 3; supp models - high AL cert).csv', row.names = F)  ## ignore the errors


WA_mean_elev = readRDS(coefs[10]) %>% pool.table(type = "all") %>% 
  mutate(mod = coefs[10], term =  factor(term, levels = c('(Intercept)', y_orders), labels = c('(Intercept)', y_labs))) %>% arrange(term)

table6 = WA_mean_elev %>% select(c(term, estimate, std.error, conf.low, conf.high, p.value, mod)) %>% num_tidy() %>% 
  write.csv('../SI Tables/Table S10 HWA results (Related to figure 3; supp models - mean elev).csv', row.names = F)  ## ignore the errors



## ............ Figure 4. forest plots by elevation band  ##############

library(mice)
library(tidyverse)

##### Part 1: Average coefficients from the 100 trees using Rubin's rule 

ave_100_by_group <- function(coef_table){
  
    coef_mean_ntree=c()
    for (bin in unique(coef_table$elevbin) ){
        for (y in unique(coef_table$Y)) {
          
            # this gives the 100 model results per elevation band per wing metric
            coef_bin_y_i <- coef_table %>% filter(elevbin == bin, Y == y)
            
            coef_bin_y_i$term = coef_bin_y_i$predictor
            coef_bin_y_i$df.residual = coef_bin_y_i$sampleSize - length(unique(coef_bin_y_i$term))
            
            colnames(coef_bin_y_i)[c(1,2)] <- c('estimate', 'std.error')
            
            pooled_bin_y_i = pool.table(coef_bin_y_i, type = "all")
            
            # add back the band info
            pooled_bin_y_i$Y = y
            pooled_bin_y_i$elevbin = bin
            pooled_bin_y_i$sampleSize = unique(coef_bin_y_i$sampleSize)
            
            coef_mean_ntree = rbind(coef_mean_ntree, pooled_bin_y_i)
        }
    }
    return(coef_mean_ntree)
}



##### Part 2: plot the slopes

plot_by_group <- function(coef){
    
  coef$sig = ifelse(coef$p.value < 0.05, 1, 0.3)
  coef$elevbin = as.character(coef$elevbin)
  
  coef_elev <<- filter(coef, term %in% c('Elev'))

  coef_elev$Y <- factor(coef_elev$Y, levels = c('log.Sec1', 'log.Wing.Length', 'log.WA', 'HWI'))

  p = ggplot(coef_elev, aes(elevbin, estimate, col = Y, alpha = I(sig))) + 
          geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                          linewidth = 1.4, 
                          position = position_dodge(width = 0.55)) +  
          geom_hline(yintercept = 0, lty = 'solid', col = 'grey80') + 
          labs(y = 'Effect size of elevation')+
          theme_classic() +
          scale_color_manual(name = '',
                    values=c('brown', 'olivedrab', '#e68932', '#192370'), 
                    guide = guide_legend(reverse = TRUE),
                    labels = c('Wing width', 'Wing length','Hand-wing area','Hand-wing index')) +
          coord_flip() +
          theme(axis.text = element_text(size=14),
                axis.title = element_text(size=16),
                legend.title = element_text(size=12),
                legend.text = element_text(size=12))
    
    return(p)
    
}


### Figure 4
  
setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/elev_band_coef100/")

pdf('../../figures/Figure 4 (4 metrics).pdf', width = 9.5, height = 5.3)
ave_100_by_group(readRDS('full (ref=flap).rds')) %>% plot_by_group() +
    scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-4\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-5\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-6\n(n=', unique(coef_elev$sampleSize)[4], ')'),
                            paste0('4-7\n(n=', unique(coef_elev$sampleSize)[5], ')'),
                            paste0('5-8\n(n=', unique(coef_elev$sampleSize)[6], ')')))
dev.off()


################# Supplementary analyses to Figure 4 ########################

## Figure S3: non-migratory landbird model #######

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/elev_band_coef100/")

pdf('../../figures/SI/Fig. S4. elev bands (land.sedentary, 4 metrics).pdf', width = 7.2, height = 5.3)
ave_100_by_group(readRDS('land.sedentary NEW (ref=flap).rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-4\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-5\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-6\n(n=', unique(coef_elev$sampleSize)[4], ')'), 
                            paste0('4-7\n(n=', unique(coef_elev$sampleSize)[5], ')'), 
                            paste0('5-8\n(n=', unique(coef_elev$sampleSize)[6], ')'))) +
  theme(legend.position = c(0.82, 0.17), legend.title = element_blank(), 
        plot.margin = margin(l=7, t=7, r=15))
dev.off()


## Figure S4: additional sensitivity analyses #########

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/elev_band_coef100/")

# using high-certainty data

pdf('../../figures/SI/Fig. S5a. elev bands (full; high.certainty, 4 metrics).pdf', width = 6.2, height = 5.2)
ave_100_by_group(readRDS('full (ref=flap) high.certainty.rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-4\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-5\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-6\n(n=', unique(coef_elev$sampleSize)[4], ')'),
                            paste0('4-7\n(n=', unique(coef_elev$sampleSize)[5], ')'),
                            paste0('5-8\n(n=', unique(coef_elev$sampleSize)[6], ')'))) +
  theme(legend.position = 'none')
dev.off()


# using mean elevation

pdf('../../figures/SI/Fig. S5b. elev bands (full; mean.elevation; 4 metrics).pdf', width = 6.2, height = 5.2)
ave_100_by_group(readRDS('full (ref=flap) mean.elevation (band division uses MEAN).rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-3\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-4\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-5\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-6\n(n=', unique(coef_elev$sampleSize)[4], ')'))) +
  theme(legend.position = c(0.82, 0.17), legend.title = element_blank(), 
        plot.margin = margin(l=7, t=7, r=15))
dev.off()



## Figure S6: using alternative bandwidth settings (2 km and 4 km) #########

# use bandwidth = 4 km

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/SI4_elev_band_coef100_BW.4km/")
pdf('../../figures/SI/Fig. S6a. bw=4km.pdf', width = 6.2, height = 5.2) 
ave_100_by_group(readRDS('full (ref=flap; bw=4km).rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-4\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-5\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-6\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-7\n(n=', unique(coef_elev$sampleSize)[4], ')'),
                            paste0('4-8\n(n=', unique(coef_elev$sampleSize)[5], ')'))) +
  theme(legend.position = 'none')
dev.off()


# use bandwidth = 2 km

setwd("C:/Users/Jinne/OneDrive - Imperial College London/_WING/elev-r/hpc coef results/SI4_elev_band_coef100_BW.2km/")
pdf('../../figures/SI/Fig. S6b. bw=2km.pdf', width = 6.2, height = 5.2)
ave_100_by_group(readRDS('full (ref=flap; bw=2km).rds')) %>% plot_by_group() +
  scale_x_discrete(name ="Elevation band (km)\n", 
                   breaks=as.character(unique(coef_elev$elevbin)),
                   labels=c(paste0('0-2\n(n=', unique(coef_elev$sampleSize)[1], ')'),
                            paste0('1-3\n(n=', unique(coef_elev$sampleSize)[2], ')'),
                            paste0('2-4\n(n=', unique(coef_elev$sampleSize)[3], ')'), 
                            paste0('3-5\n(n=', unique(coef_elev$sampleSize)[4], ')'),
                            paste0('4-6\n(n=', unique(coef_elev$sampleSize)[5], ')'),
                            paste0('5-7\n(n=', unique(coef_elev$sampleSize)[6], ')'))) +
  theme(legend.position = c(0.82, 0.17), legend.title = element_blank(), 
        plot.margin = margin(l=7, t=7, r=15))
dev.off()



