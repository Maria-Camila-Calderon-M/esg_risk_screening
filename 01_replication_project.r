## University of Geneva
## Course: Environmental economics and climate change
## Fall semester 2025
## Group 5
## Andrea Acosta, Catalina Munar, Mara Moos, Jennifer Yañez, María Camila Calderón
## This code is divided in two parts: 
## Part 1 - Replication
## Part 2 - Extension
## Last updated: January 6th 2026

#### REPLICATION ####
## This code conducts analysis for reforestation 
## dynamics in indigenous lands in Brazil
## Authors: Kathryn Baragwanath, Ella Bayi and Nilesh Shinde
## Created:  July 2022
## Last Updated: 20 April 2023

# Begin this project in a fresh environment
install.packages("purrr")

#package ‘hrbr’ is not available for the team's version of R
#package ‘rgdal’ is not available for the team's version of R

# Loading Packages
library(tidyverse)
library(sf)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
require(maps)
library(readxl)
library(stringr)
library(scales)
library(readr)
library(data.table)
library(lattice)
library(viridis)
library(ggthemes)
library(rasterVis)
library(RColorBrewer)
library(rdrobust)
library(fixest)
library(did)
library(ggpubr)
library(grid)
library(ggplotify)
library(ggplot2)
library(ggpubr)
library(modelsummary)
library(here)

## Setting seed for replication
set.seed(1234)

## Loading data
load("points_final.Rdata")

## Modify code for tables to Latex
tidy.rdrobust <- function(x, ...) {
  ret <- data.frame(
    term = row.names(x$coef),
    estimate = x$coef[, 1],
    std.error = x$se[, 1],
    p.value = x$pv[, 1]
  )
  row.names(ret) <- NULL
  ret
}

glance.rdrobust <- function(x, ...) {
  ret <- data.frame(
    "Mean Control" = x$beta_p_l[1],
    Kernel = x$kernel,
    Bandwidth = x$bwselect,
    N = as.integer(round(x$N_h[1] + x$N_h[2], digits=1)),
    BW = x$bws[1,1]
  )
  ret
}

####################################################
####################################################
##############         RDD          ################
####################################################
####################################################

### 15 years After Homologation 

p <- points %>%   group_by(nearit_FID,year) %>%
  mutate(t_id = max(t)) %>% subset(t_id==1 & potential_ref>=0.027) 

p<-subset(p, year>=Date & (year<= Date +15))

p<-subset(p,  (near_line >= -20000 & near_line <= 20000))
p<-subset(p, near_line<=-1000 | near_line >=1000)
p$near_line<- ifelse(p$near_line<0, p$near_line+999,p$near_line-999)

#check the data#
View(p)

Y <- p$extent_prop
X <- p$near_line

rdplot_extent_after <-rdplot(Y,X,c=0,  p=1,  nbins=100, h=12000, x.lim = c(-12000,12000), 
                             kernel="tri", y.lim = c(0,35), title = "" ,x.label = "Distance to Border (m)",
                             y.label = "Percent Extent")

rdplot_extent_after<- rdplot_extent_after$rdplot + theme(axis.text=element_text(size=16),
                            axis.title=element_text(size=16))

png("rd_extent.png", width = 4, height = 4, units='in', res=300)
rdplot_extent_after
dev.off()

mod1_extent_after<-rdrobust(Y,X, c=0, p=1, kernel="triangular",vce="nn", cluster = p$nearit_FID)
summary(mod1_extent_after)

## Average age of Secondary Forest 

Y <- p$mean_age
X <- p$near_line

rdplot_age_after <-rdplot(Y,X,c=0,  p=1,  nbins=100, h=12000, x.lim = c(-12000,12000), y.lim = c(0,15) , kernel="tri",  title = "" ,x.label = "Distance to Border (m)", y.label = "Average Age of Secondary Forest")

rdplot_age_after<- rdplot_age_after$rdplot + theme(axis.text=element_text(size=16),
                           axis.title=element_text(size=16))

png("rd_age.png", width = 4, height = 4, units='in', res=300)
rdplot_age_after
dev.off()

mod1_age_after<-rdrobust(Y,X, c=0, p=1, kernel="triangular",vce="nn", cluster = p$nearit_FID)
summary(mod1_age_after)

### 15 years Before Homologation 

p<- subset(points, !is.na(bef_hom_fid))
p<-subset(p, (year<Date & year>=Date-15)  & potential_ref>=0.027)
p<-subset(p,  (near_line >= -20000 & near_line <= 20000))
p<-subset(p, near_line<=-1000 | near_line >=1000)
p$near_line<- ifelse(p$near_line<0, p$near_line+999,p$near_line-999)

#RDD for Silva Junior Secondary Forest Extent
Y <- p$extent_prop
X <- p$near_line

rdplot_extent_bef <-rdplot(Y,X,c=0, nbins=100, p=1, h=12000, x.lim = c(-12000,12000), kernel="tri", y.lim = c(0,35), title = "" ,x.label = "Distance to Border (m)", y.label = "Percent Extent")

rdplot_extent_bef <- rdplot_extent_bef$rdplot + theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

png("rd_extent.png", width = 4, height = 4, units='in', res=300)
rdplot_extent_bef
dev.off()

mod1_extent_bef<-rdrobust(Y,X, c=0, p=1, kernel="triangular", vce = "nn", cluster = p$nearit_FID)
summary(mod1_extent_bef)

### Average Age of Secondary Forest 

Y <- p$mean_age
X <- p$near_line

rdplot_age_bef <-rdplot(Y,X,c=0, nbins=100, p=1, h=12000, x.lim = c(-12000,12000), y.lim = c(0,15), kernel="tri",  title = "" ,x.label = "Distance to Border (m)", y.label = "Average Age of Secondary Forest")

rdplot_age_bef <- rdplot_age_bef$rdplot + theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

png("rd_age.png", width = 4, height = 4, units='in', res=300)
rdplot_age_bef
dev.off()

mod1_age_bef<-rdrobust(Y,X, c=0, p=1, kernel="triangular", vce = "nn", cluster = p$nearit_FID)
summary(mod1_age_bef)

### Non Homologated

p<-subset(points, !is.na(nhom_fid) & potential_ref>=0.027)
p<-subset(p,  (near_line >= -20000 & near_line <= 20000) )
p<-subset(p, near_line<=-1000 | near_line >=1000)
p$near_line<- ifelse(p$near_line<0, p$near_line+999,p$near_line-999)

#RDD for Silva Junior Secondary Forest Extent

Y <- p$extent_prop
X <- p$near_line

rdplot_extent_nhom <-rdplot(Y,X,c=0, nbins=150, p=1, h=12000, x.lim = c(-12000,12000), kernel="tri", y.lim = c(0,35), title = "" ,x.label = "Distance to Border (m)", y.label = "Percent Extent")

rdplot_extent_nhom <- rdplot_extent_nhom$rdplot + theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

png("rd_extent.png", width = 4, height = 4, units='in', res=300)
rdplot_extent_nhom
dev.off()

mod1_extent_nhom<-rdrobust(Y,X, c=0, p=1, kernel="triangular", vce = "nn")
summary(mod1_extent_nhom)

# Average Age of Secondary Forest

Y <- p$mean_age
X <- p$near_line

rdplot_age_nhom <-rdplot(Y,X,c=0, nbins=100, p=1, h=12000, x.lim = c(-12000,12000), y.lim=c(0,15) , kernel="tri",  title = "" ,x.label = "Distance to Border (m)", y.label = "Average Age of Secondary Forest")

rdplot_age_nhom <- rdplot_age_nhom$rdplot + theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

png("rd_age.png", width = 4, height = 4, units='in', res=300)
rdplot_age_nhom
dev.off()

mod1_age_nhom<-rdrobust(Y,X, c=0, p=1, kernel="triangular", vce = "nn")
summary(mod1_age_nhom)

## Rdplots 
# Comparing Rdd Plots code

a<-rdplot_extent_nhom
b<-rdplot_extent_bef
c<-rdplot_extent_after

d<-rdplot_age_nhom
e<-rdplot_age_bef
f<-rdplot_age_after


plot<-ggarrange(
  a,b,c,
  nrow = 1, ncol = 3,
  labels = c("Non Hom", "Before", "After")
)
png("rd_comparison_extent_sj.png", width = 15, height = 6, units='in', res=300 )
plot
dev.off()

plot<-ggarrange(
  d,e,f,
  nrow = 1, ncol = 3,
  labels = c("Non Hom", "Before", "After")
)
png("rd_comparison_age_sj.png", width = 15, height = 6, units='in', res=300 )
plot
dev.off()

##### Tables #####
## Extent

models<-list(
  "Non Homologated"= mod1_extent_nhom,
  "Before Homologation" = mod1_extent_bef,
  "After Homologation" = mod1_extent_after
)

modelsummary(models, estimate = "{estimate}{stars}", output = "html", title = "Secondary Vegetation, Silva Junior")
#We don't use Latex, we use html and viewer

## From this output we keep only the robust estimates. 

## AGE 

models<-list(
  "Non Homologated"= mod1_age_nhom,
  "Before Homologation" = mod1_age_bef,
  "After Homologation" = mod1_age_after
)

modelsummary(models, estimate = "{estimate}{stars}", output = "html", title = "Secondary Vegetation Age, Silva Junior")

## From this output we keep only the robust estimates. 


########  COEFPLOTS FOR DIFFERENT OUTCOME VARIABLES #########

###### Extent SJ ######
#After Homologation 

coef<-mod1_extent_after$coef
se<-mod1_extent_after$se

dat1 <- data.frame(coef, se)
dat1$type <- c("After Homologation", "After Homologation", "After Homologation")
dat1$type2<- c("Conventional", "Bias", "Robust")
dat1$val<-2


# Before Homologation 
coef<-mod1_extent_bef$coef
se<-mod1_extent_bef$se

dat2 <- data.frame(coef, se)
dat2$type <- c("Before Homologation", "Before Homologation", "Before Homologation")
dat2$type2<- c("Conventional", "Bias", "Robust")
dat2$val<-4


#Non Homologated
coef<-mod1_extent_nhom$coef
se<-mod1_extent_nhom$se

dat3 <- data.frame(coef, se)
dat3$type <- c("Non Homologated", "Non Homologated", "Non Homologated")
dat3$type2<- c("Conventional", "Bias", "Robust")
dat3$val<-6


dat<-rbind(dat1, dat2,  dat3)

dat<-subset(dat, type2=="Robust")

plot<-data.frame(dat,stringsAsFactors = F)
names(plot)<-c("Estimate","S.E.","type","test", "val")
plot$Estimate<-as.numeric(as.character(plot$Estimate))
plot$S.E.<-as.numeric(as.character(plot$S.E.))

#get upper and lower CI bounds
plot$lower<-plot$Estimate-1.96*plot$S.E.
plot$upper<-plot$Estimate+1.96*plot$S.E.
plot$lower1<-plot$Estimate-1.645*plot$S.E.
plot$upper1<-plot$Estimate+1.645*plot$S.E.
type<-plot$type
plot$type<-factor(plot$type,levels = plot$type[order(plot$val)])

#Adjust sizes and transparencies
plot$main<-c(rep("D",3))
head(plot)

mytheme <- theme(
  axis.line = element_line(linewidth = 1, colour = "black"),
  panel.background = element_rect(fill = "white")
)

plot1 <- ggplot() +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype=1, color="black", linewidth=0.8) +
  
  # Add point estimates
  geom_point(aes(x = Estimate, y = type), color="grey", size=2, data = plot) +
  
  # Confidence intervals
  geom_segment(aes(y = type, x = lower, yend = type, xend = upper), color="grey", linewidth=0.8, data = plot) +
  
  # Text and theme
  theme(axis.text.y = element_text(size = 18, angle = 0, hjust = 1),
        axis.text.x = element_text(size = 18),
        text = element_text(size = 18, colour = "black"),
        legend.position = "none") +
  
  # Color scales and custom theme
  scale_color_manual(values = c("black", "grey")) +
  mytheme +
  
  # Labels
  xlab("SV Extent") +
  ylab("")

# Display the plot Extent
plot1

## We get the second plot 

###### Age SJ ######
#After Homologation 

coef<-mod1_age_after$coef
se<-mod1_age_after$se

dat1 <- data.frame(coef, se)
dat1$type <- c("After Homologation", "After Homologation", "After Homologation")
dat1$type2<- c("Conventional", "Bias", "Robust")
dat1$val<-2

# Before Homologation 
coef<-mod1_age_bef$coef
se<-mod1_age_bef$se

dat2 <- data.frame(coef, se)
dat2$type <- c("Before Homologation", "Before Homologation", "Before Homologation")
dat2$type2<- c("Conventional", "Bias", "Robust")
dat2$val<-4

#Non Homologated
coef<-mod1_age_nhom$coef
se<-mod1_age_nhom$se

dat3 <- data.frame(coef, se)
dat3$type <- c("Non Homologated", "Non Homologated", "Non Homologated")
dat3$type2<- c("Conventional", "Bias", "Robust")
dat3$val<-6


dat<-rbind(dat1, dat2,  dat3)

dat<-subset(dat, type2=="Robust")

plot<-data.frame(dat,stringsAsFactors = F)
names(plot)<-c("Estimate","S.E.","type","test", "val")
plot$Estimate<-as.numeric(as.character(plot$Estimate))
plot$S.E.<-as.numeric(as.character(plot$S.E.))

#get upper and lower CI bounds
plot$lower<-plot$Estimate-1.96*plot$S.E.
plot$upper<-plot$Estimate+1.96*plot$S.E.
plot$lower1<-plot$Estimate-1.645*plot$S.E.
plot$upper1<-plot$Estimate+1.645*plot$S.E.
type<-plot$type
plot$type<-factor(plot$type,levels = plot$type[order(plot$val)])
#Adjust sizes and transparencies
plot$main<-c(rep("D",3))
head(plot)

# Corrected theme
mytheme <- theme(
  axis.line = element_line(linewidth = 1, colour = "black"),
  panel.background = element_rect(fill = "white")
)

# Create plot2
plot2 <- ggplot() +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype = 1, color = "black", linewidth = 0.8) +
  
  # Add point estimates
  geom_point(aes(x = Estimate, y = type), color = "grey", size = 2, data = plot) +
  
  # Confidence intervals
  geom_segment(aes(y = type, x = lower, yend = type, xend = upper), color = "grey", linewidth = 0.8, data = plot) +
  
  # Text and theme
  theme(
    axis.text.y = element_text(size = 18, angle = 0, hjust = 1),
    axis.text.x = element_text(size = 18),
    text = element_text(size = 18, colour = "black"),
    legend.position = "none"
  ) +
  
  # Optional color scale and custom theme
  scale_color_manual(values = c("black", "grey")) +
  mytheme +
  
  # Labels
  xlab("SV Age") +
  ylab("")

# Display the plot Age
plot2

# Now we merge both and download as png

plot<-ggarrange(
  plot1, plot2,
  # labels = c("SV Extent (MB)", "SV Extent (SJ)", "SV Increment (MB)", "SV Increment (SJ)"),
  nrow = 1, ncol = 2
)

png("coefplot_sv_age.png",width = 12, height = 4, units='in', res=300 )
print(plot)
dev.off()


####################################################
####################################################
##############   EVENT STUDY TWFE   ################
####################################################
####################################################

## Including model structure into modelsummary

tidy.AGGTEobj<- function(x, ...){
  ret <- data.frame(
    term = c("ATT"),
    estimate = x[["overall.att"]],
    std.error = x[["overall.se"]]
  )
  ret
}

glance.AGGTEobj <- function(x, ...) {
  ret <- data.frame(
    "Model" = x[["type"]],
    "nobs" = x[["DIDparams"]][["n"]],
    "ntime" =x[["DIDparams"]][["nT"]],
    "vcov.type" = x[["DIDparams"]][["clustervars"]]
    
  )
  ret
}


## Getting data ready 

p<-subset(points,  (near_line >= -20000 & near_line <= 20000))
p<-subset(p, near_line<=-1000 | near_line >=1000)
p$near_line<- ifelse(p$near_line<0, p$near_line+999,p$near_line-999)

library(data.table)
p<-as.data.table(p)
p[, time_to_treat := ifelse(treatment==1, year - `Date`, NA)]

library(tidyr)

p <- p %>%   group_by(nearit_FID) %>%
  fill(time_to_treat) %>%
  fill(time_to_treat, .direction = "up")

## Dropping always treated units from sample. 
p$g<-ifelse(p$treatment==1, p$Date, 0)
table(p$g)

p_sub<-subset(p,g==0 | g>=1989)

library(did)

did.att.gt.cov <- att_gt(yname="extent_prop",
                        tname="year",
                        idnam="FID_1",
                        gname="g",
                        control_group = "notyettreated",
                        clustervars = "nearit_FID",
                        xformla = ~ near_mines + near_roads + elevation1,
                        data=p_sub)

# summary(did.att.gt)

es_dynamic<-aggte(did.att.gt.cov, type = "dynamic")
summary(es_dynamic)

dynamic.plot.cov<- ggdid(es_dynamic, xgap = 5, title=" ", ) +
  theme(axis.text.y = element_text(size = 18, angle = 0, hjust = 1),
        axis.text.x = element_text(size = 18)) + theme(legend.text=element_text(size=18)) +
  xlim(-10,25) + ylim(-1.5,9) 

dynamic.plot.cov

png("event_did_dynamic.png",width = 10, height = 6, units='in', res=300 )
dynamic.plot.cov
dev.off()

es_simple<-aggte(did.att.gt.cov, type = "simple")
summary(es_simple)

es_group <- aggte(did.att.gt.cov, type = "group")
summary(es_group)

es_calendar <- aggte(did.att.gt.cov, type = "calendar")
summary(es_calendar)

sv_simple<-data.frame(es_simple[["overall.att"]], es_simple[["overall.se"]], es_simple[["type"]],  es_simple[["DIDparams"]][["n"]], es_simple[["DIDparams"]][["nT"]])
colnames(sv_simple)<-c("coeff", "std.error", "type", "n", "periods")

sv_dynamic<-data.frame(es_dynamic[["overall.att"]], es_dynamic[["overall.se"]], es_dynamic[["type"]], es_dynamic[["DIDparams"]][["n"]], es_dynamic[["DIDparams"]][["nT"]])
colnames(sv_dynamic)<-c("coeff", "std.error", "type", "n", "periods")

sv_calendar<-data.frame(es_calendar[["overall.att"]], es_calendar[["overall.se"]], es_calendar[["type"]], es_calendar[["DIDparams"]][["n"]], es_calendar[["DIDparams"]][["nT"]])
colnames(sv_calendar)<-c("coeff", "std.error", "type", "n", "periods")

sv_group<-data.frame(es_group[["overall.att"]], es_group[["overall.se"]], es_group[["type"]], es_group[["DIDparams"]][["n"]], es_group[["DIDparams"]][["nT"]])
colnames(sv_group)<-c("coeff", "std.error", "type", "n", "periods")

sv_res_extent<-rbind( sv_dynamic,sv_simple, sv_calendar, sv_group)
sv_res_extent

### Age of secondary forest
did.att.gt.age <- att_gt(yname="mean_age",
                         tname="year",
                         idnam="FID_1",
                         gname="g",
                         clustervars = "nearit_FID",
                         control_group = "notyettreated",                      
                         data=p_sub)

# summary(did.att.gt)

es_dynamic.cov.age<-aggte(did.att.gt.age, type = "dynamic")
summary(es_dynamic.cov.age)

dynamic.plot.cov.age<- ggdid(es_dynamic.cov.age, xgap = 5, title=" ", ) +
  theme(axis.text.y = element_text(size = 18, angle = 0, hjust = 1), 
        axis.text.x = element_text(size = 18)) + theme(legend.text=element_text(size=18)) +
  xlim(-10,25) + ylim(-1,12)

dynamic.plot.cov.age

png("event_did_dynamic_age.png",width = 10, height = 6, units='in', res=300 )
dynamic.plot.cov.age
dev.off()

## Table 
es_simple_age<-aggte(did.att.gt.age, type = "simple")
summary(es_simple_age)

es_group_age <- aggte(did.att.gt.age, type = "group")
summary(es_group_age)

es_calendar_age <- aggte(did.att.gt.age, type = "calendar")
summary(es_calendar_age)


sv_simple<-data.frame(es_simple_age[["overall.att"]], es_simple_age[["overall.se"]], es_simple_age[["type"]],  es_simple_age[["DIDparams"]][["n"]], es_simple_age[["DIDparams"]][["nT"]])
colnames(sv_simple)<-c("coeff", "std.error", "type", "n", "periods")

sv_dynamic<-data.frame(es_dynamic.cov.age[["overall.att"]], es_dynamic.cov.age[["overall.se"]], es_dynamic.cov.age[["type"]], es_dynamic.cov.age[["DIDparams"]][["n"]], es_dynamic.cov.age[["DIDparams"]][["nT"]])
colnames(sv_dynamic)<-c("coeff", "std.error", "type", "n", "periods")

sv_calendar<-data.frame(es_calendar_age[["overall.att"]], es_calendar_age[["overall.se"]], es_calendar_age[["type"]], es_calendar_age[["DIDparams"]][["n"]], es_calendar_age[["DIDparams"]][["nT"]])
colnames(sv_calendar)<-c("coeff", "std.error", "type", "n", "periods")

sv_group<-data.frame(es_group_age[["overall.att"]], es_group_age[["overall.se"]], es_group_age[["type"]], es_group_age[["DIDparams"]][["n"]], es_group_age[["DIDparams"]][["nT"]])
colnames(sv_group)<-c("coeff", "std.error", "type", "n", "periods")

sv_res<-rbind( sv_dynamic, sv_simple, sv_calendar, sv_group)

sv_res

plot<-ggarrange(
  dynamic.plot.cov , dynamic.plot.cov.age,
  labels = c("A. SV Extent", "B. SV Age"),
  nrow = 1, ncol = 2
)

png("event_did.png",width = 12, height = 6, units='in', res=300 )
plot
dev.off()

#### EXTENSION ####

## TESTING H1:Secondary forest growth is bigger with an increasing distance to mines

#Takaki's notes about the points
p %>% count(FID_1)
p %>% distinct(FID_1)
p %>% ungroup() %>%  distinct(FID_1)

#### Option 1: Linear effect ####
## Linear effect of distance to mines to secondary forest extent and age

# With controls roads and elevation
mod_extent_main <- feols(extent_prop ~ near_mines + near_roads + elevation1 | FID_1 + year, data = points) 
mod_age_main  <- feols(mean_age ~ near_mines + near_roads + elevation1 | FID_1 + year, data = points) 
#The model automatically removes roads and elevation because of collinearity

# Robustness checks (one control at a time) 
# Only roads
mod_extent_roads <- feols(extent_prop ~ near_mines + near_roads | FID_1 + year, 
                          data = points, cluster = ~FID_1)
#The model automatically removes roads because of collinearity

# Only elevation
mod_extent_elev  <- feols(extent_prop ~ near_mines + elevation1 | FID_1 + year, 
                          data = points, cluster = ~FID_1)
#The model automatically removes elevation because of collinearity

#Tables
models_extent <- list(
  "Linear FE" = mod_extent_main,
  "Robustness: + Roads" = mod_extent_roads,
  "Robustness: + Elevation" = mod_extent_elev
)

models_age <- list(
  "Main: Linear FE" = mod_age_main)

#Extent

modelsummary(models_extent,
             stars = TRUE,                    # significance stars
             statistic = "std.error",         # show SEs below coefficients
             coef_rename = c("near_mines" = "Distance to Mines (km)"),
             gof_omit = 'Adj|Within',         # hide Adj/Within R² if desired
             output = "linear_FE_extent.html")  # saves as HTML

#Age
modelsummary(models_age,
             stars = TRUE,
             statistic = "std.error",
             coef_rename = c("near_mines" = "Distance to Mines (km)"),
             gof_omit = 'Adj|Within',
             output = "linear_FE_age.html")

models <- list(
  "Secondary Forest Extent" = mod_extent_main,
  "Secondary Forest Age"   = mod_age_main
)

modelsummary(models,
             stars = TRUE,                    # significance stars
             statistic = "std.error",         # show standard errors below coefficients
             coef_rename = c("near_mines" = "Distance to Mines (km)"),
             gof_omit = 'Adj|Within',         # omit Adj/Within R² if desired
             output = "linear_FE_extent_age.html")

#### Option 2: No linear - Bins ####

#Explore near_mines

# Summary and histogram
summary(points$near_mines)

hist1 <- ggplot(points, aes(x = near_mines)) +
  geom_histogram(binwidth = 10000, fill = "steelblue", color = "white") +
  scale_y_continuous(labels = comma) +  # plain numbers with commas
  labs(title = "Distribution of distance to mines (km)",
       x = "Distance to mines (km)",
       y = "Count")

ggsave("histogram_near_mines.png", plot = hist1, width = 8, height = 5, dpi = 300)

library(scales)

# Histogram with vertical lines for quartiles
hist2 <-ggplot(points, aes(x = near_mines)) +
  geom_histogram(binwidth = 10000, fill = "steelblue", color = "white") +
  scale_y_continuous(labels = comma) +  # show plain numbers
  labs(title = "Distribution of distance to mines (km)",
       x = "Distance to mines (km)",
       y = "Count") +
  geom_vline(xintercept = c(62502.5, 118553.1, 215192.7), 
             color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = c(62502.5, 118553.1, 215192.7),
           y = max(hist(points$near_mines, plot = FALSE)$counts)*0.95,
           label = c("Q1", "Median", "Q3"),
           angle = 90, vjust = -0.5, color = "red", size = 3)

ggsave("histogram_near_mines_quartiles.png", plot = hist2, width = 8, height = 5, dpi = 300)

##Option 2.1. Bins according to histogram
# Create bins
points$near_mines_bin <- cut(points$near_mines,
                             breaks = c(0, 50000, 120000, 215000, 400000, Inf),
                             labels = c("Very Close", "Close", "Moderate", "Far", "Very Far"),
                             include.lowest = TRUE,
                             right = TRUE)

# Check distribution of observations in bins
table(points$near_mines_bin)

library(fixest)

mod_extent_bins <- feols(
  extent_prop ~ near_mines_bin | year,
  data = points,
  cluster = ~FID_1
)

mod_age_bins <- feols(
  mean_age ~ near_mines_bin | year,
  data = points,
  cluster = ~FID_1
)

library(modelsummary)

models_bins <- list(
  "Secondary Forest Extent" = mod_extent_bins,
  "Secondary Forest Age"   = mod_age_bins
)

modelsummary(
  models_bins,
  stars = TRUE,
  statistic = "std.error",
  coef_rename = c(
    "near_mines_binVery Close" = "Very Close (0–50k km)",
    "near_mines_binClose"      = "Close (50k–120k km)",
    "near_mines_binModerate"   = "Moderate (120k–215k km)",
    "near_mines_binFar"        = "Far (215k–400k km)",
    "near_mines_binVery Far"   = "Very Far from mines (>400 km)" 
  ),
  gof_map = c(
    "nobs", "r.squared", "rmse", "aic", "bic"
  ),
  notes = c(
    "Reference category: Very Close (0-50k km)",
    "All models include year fixed effects.",
    "Standard errors clustered at the territory (FID_1) level."
  ),
  output = "bins_distance_to_mines_extent_age_1.html"
)

#Plots
#Coefficients with confidence intervals

tidy_extent <- broom::tidy(mod_extent_bins, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(outcome = "Secondary Forest Extent")

tidy_age <- broom::tidy(mod_age_bins, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(outcome = "Secondary Forest Age")

# Combine
tidy_combined <- bind_rows(tidy_extent, tidy_age)

# Rename bins for clarity
tidy_combined <- tidy_combined %>%
  mutate(term = recode(term,
                       "near_mines_binClose" = "Close (50k–120k km)",
                       "near_mines_binModerate" = "Moderate (120k–215k km)",
                       "near_mines_binFar" = "Far (215k–400k km)",
                       "near_mines_binVery Far" = "Very Far (>400 km)"),
         term = factor(term, levels = c("Close (50k–120k km)",
                                        "Moderate (120k–215k km)",
                                        "Far (215k–400k km)",
                                        "Very Far (>400 km)")))

# Plot
p <- ggplot(tidy_combined, aes(x = estimate, y = term, color = outcome)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Estimated Effect", y = "Distance to Mines (bins)", color = "") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("steelblue", "darkorange")) +
  theme(legend.position = "top")

# Show plot
print(p)

# Export as PNG
ggsave("bins_coefficients_combined.png", plot = p,
       width = 10, height = 5, dpi = 300)

##Option 2.2. Bins according to quantiles

#Plot actual distributions by bin
plot1 <- ggplot(points, aes(x = near_mines_qbin, y = extent_prop)) +
  geom_boxplot() +
  labs(title = "Secondary Forest Extent by Distance Quantile Bin")

ggsave(
  filename = "secondary_forest_extent_by_bin.png",  # name of the file
  plot = plot1,                                        # plot object
  width = 8, height = 6,                           # size in inches
  dpi = 300                                        # resolution
)

#Create bins
library(dplyr)

#Create quantiles (1 = closest, 5 = farthest)
points$near_mines_qbin <- cut(
  points$near_mines,
  breaks = quantile(points$near_mines, probs = seq(0, 1, 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Bin 1 (Closest)", "Bin 2", "Bin 3", "Bin 4", "Bin 5 (Farthest)")  # 5 bins
)

#Check quantiles
table(points$near_mines_qbin)
quantile_breaks <- quantile(points$near_mines, probs = seq(0, 1, 0.2), na.rm = TRUE)
quantile_breaks

#Check labels
breaks <- c(110.625, 52461.450, 94104.280, 148706.900, 246653.400, 1005432.326)

breaks_formatted <- c(
  "111 km",        # round 110.625
  "52k km",        # 52461.450
  "94k km",        # 94104.280
  "149k km",       # 148706.900
  "247k km",       # 246653.400
  "1005k km"       # 1005432.326
)

# Create bin labels
bin_labels <- paste0(breaks_formatted[-length(breaks_formatted)], " – ",
                     breaks_formatted[-1])

bin_labels

points$near_mines_qbin <- cut(
  points$near_mines,
  breaks = breaks,
  include.lowest = TRUE,
  labels = bin_labels
)

table(points$near_mines_qbin)

#Let's try the models with three different bins as reference

run_models_qbin <- function(ref_bin_label) {
  
  # Relevel factor
  points$near_mines_qbin <- relevel(points$near_mines_qbin, ref = ref_bin_label)
  
  # Secondary Forest Extent
  mod_extent <- feols(
    extent_prop ~ near_mines_qbin | factor(year),
    data = points,
    cluster = ~FID_1
  )
  
  # Secondary Forest Age
  mod_age <- feols(
    mean_age ~ near_mines_qbin | factor(year),
    data = points,
    cluster = ~FID_1
  )
  
  list(
    "Secondary Forest Extent" = mod_extent,
    "Secondary Forest Age" = mod_age
  )
}

# Farthest bin reference
models_ref_farthest <- run_models_qbin("247k km – 1005k km")

# Closest bin reference
models_ref_closest <- run_models_qbin("111 km – 52k km")

# Median bin reference
models_ref_median <- run_models_qbin("94k km – 149k km")

# Combine all models in one list
quantile_models <- c(
  models_ref_farthest,
  models_ref_closest,
  models_ref_median
)

# Optional: create nicer names for the columns
model_names <- c(
  "Extent (Farthest Ref)", "Age (Farthest Ref)",
  "Extent (Closest Ref)", "Age (Closest Ref)",
  "Extent (Median Ref)",  "Age (Median Ref)"
)

# Export table
modelsummary(
  quantile_models,
  stars = TRUE,
  statistic = "std.error",
  coef_rename = c(
    "near_mines_qbin111 km – 52k km" = "111 km – 52k km (Closest)",
    "near_mines_qbin52k km – 94k km" = "52k km – 94k km",
    "near_mines_qbin94k km – 149k km" = "94k km – 149k km (Median)",
    "near_mines_qbin149k km – 247k km" = "149k km – 247k km",
    "near_mines_qbin247k km – 1005k km" = "247k km – 1005k km (Farthest)"
  ),
  gof_map = c("nobs", "r.squared", "rmse"),
  notes = c(
    "Models run with different reference bins: farthest, closest, median.",
    "All models include year fixed effects.",
    "Standard errors clustered at the territory (FID_1) level."
  ),
  output = "bins_quantile_distance_to_mines_combined.html",
  model_names = model_names
)

####################################################################
########################## The end #################################
####################################################################