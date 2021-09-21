# 25/5/2015

##### The New Statistics with R, Chapter 10: Linear Mixed-effects Models

# Note (5-2-2015) lme4 updated since script was last edited

### Metadata
 # Source: Hector et al. 1999 Science ; Spehn et al. 2005 Ecological Monographs
 # www.esapubs.org/archive/mono/M075/001
 # Aim: Effects of diversity on yield within 8 European grasslands
 # 'data.frame':	480 obs. of  9 variables:
 # Plot      : Factor w/ 480 levels (field plots = experimental units)
 # Site      : Factor w/ 8 levels "Germany","Greece" (Fieldsites; Random grouping factor)
 # Block     : Factor w/ 15 levels (2 blocks per site except portugal; Random grouping factor)
 # Mix       : Factor w/ 200 levels (partially crossed species mixtures within sites)
 # Mix_Nested: Factor w/ 220 levels (assuming full nesting of ecotypes within sites)
 # Diversity : No. of plant species per plot (primary explanatory variable (fixed effect))
 # Shoot2    : Aboveground biomass harvested in year 2 (g m-2) [Response]
 # Shoot3    : Aboveground biomass harvested in year 3 (g m-2) [Response]
 # Diversity2: Diversity log transformed to base 2

### Set up
# Clear workspace
rm(list=ls(all=TRUE)) 
# Turn off significance stars
options(show.signif.stars= FALSE) 
# set working directory (customise to your own settings)
# setwd(...)

## Packages
library(ggplot2)
# library(lme4) 
# library(arm) 

############
# Box 10.1 #
############

# Load data from file NewStats_Chapter10_Data_Biodepth.txt"
biodepth <- read.table(file.choose(), header= TRUE) 

# Convert mixtures to factor
biodepth$Mix <- factor(biodepth$Mix) 

# transform log base 2
biodepth$Diversity2 <- log(biodepth$Diversity,2) 

# Box 10.2:
str(biodepth) 

# Take year 2 as response variable
biodepth$Yield <- biodepth$Shoot2 

# Fig 10.1:
library(ggplot2)
# qplot of Yield vs number of species for a 'no pooling analysis'
qplot(Diversity, Yield, data= biodepth, geom= c("point", "smooth"), method= "lm",
      xlab= "Number of species", ylab= "Yield", main= "No pooling analysis",
      facets=.~Site)+scale_x_continuous(trans= "log2")+theme_bw()
# ggsave(file= "Fig10-1.tiff")
dev.off()

# Fig 10.2:
# qplot of Yield vs number of species for a 'complete pooling analysis'
qplot(Diversity, Yield, data= biodepth, 
      geom= c("point", "smooth"), method= "lm", main= "Complete pooling analysis",
      xlab= "Number of species", ylab= "Yield")+scale_x_continuous(trans= "log2")+theme_bw()
# ggsave(file= "Fig10-2.tiff")
dev.off()

# An incomplete Least Squares ANOVA (not fitted)
# ls0 <- lm( terms(Yield~Site+Diversity2+Site:Diversity2, keep.order= TRUE), data= biodepth)

#################### Fig. 10.3: Supplementary R code ############################
## Including species composition
# A new dataframe of two blocks at one (Swiss) site: 
CH <- subset(biodepth, Site=="Switzerland")
# Reorder by levels of diversity
CH[order( factor(CH$Diversity) ), ]
# create colour and shape combinations to specify 32 mixtures
cols <- rep(c('black','red','green4','blue'), each = 8)
shapes <- rep(c(1:8),4)
FigTen3 <- qplot(Diversity, Yield, data= CH, geom= c("point"),
                 xlab= "Number of species", ylab= "Yield", main= "",
                 facets=.~Block, colour=Mix, shape=Mix)+
                 scale_x_continuous(trans= "log2")+theme_bw()
FigTen3 + theme(legend.position="none") + scale_colour_manual(values=cols) + scale_shape_manual(values=shapes)
# ggsave(file= "Fig10-3.tiff")
dev.off()
#########

# Least Squares 'mixed model ANOVA'
# Use terms function to keep sequential order as in the model formula

ls1 <- lm( terms(Yield~Site+Block+Diversity2+Site:Diversity2+Mix+Site:Mix, keep.order= TRUE), data= biodepth)

# par(mfrow= c(2,2) )
# plot(ls1) # increasing variance and non-normality
# library(MASS)
# boxcox(ls1) # suggests square-root transformation
# But biological justification for this transformation not clear...
# ...and results are qualitatively the same.

summary(ls1) # Table of coefficients is monstrous and hard to interpret!

anova(ls1) # Warning: Ignore F and P values - some tests are vs wrong error terms!

############
# Box 10.3 #
############

### Linear mixed-effects model analysis using lme function

# nlme package:
library(nlme)

# random intercepts model
me1 <- lme(fixed= Yield~Diversity2, random= ~1|Site/Block/Mix, na.action= na.exclude, data= biodepth) 

summary(me1) 
# grouping structure of data not understood for mix:

# Number of Observations: 464
# Number of Groups: 
#                     Site          Block %in% Site Mix %in% Block %in% Site 
#                       8                       15                      429 

## Extractor functions: 
# anova(me1)     # 
# fixef(me1)     # fixed intercept and slope
# coef(me1)      # predicted intercepts and sr.log2 slopes per site 
# ranef(me1)     # ...relative to fixed intercepts
# intervals(me1) # 95% CIs for fixed effects and variance components

# detach nlme to avoid conflicts with lme4 package
detach(package: nlme) 
     
library(lme4) # lme4 package (handles crossed random effects)

# Individual site regressions using lmList:
lmList(Yield~Diversity2|Site, data= biodepth)

# Model building random effects structure using likelihood ratio tests

# Random (Diversity) intercepts AND slopes for Sites
mer1 <- lmer(Yield~Diversity2 +(1+Diversity2|Site)+(1|Block)+(1|Mix)+(1|Site:Mix), data= biodepth) 

# check grouping has been properly understood
summary(mer1) 

# Test site:mix interaction:
mer2 <- lmer(Yield~Diversity2 +(1+Diversity2|Site)+(1|Block)+(1|Mix)             , data= biodepth) 

# more complex mer0 is preferred - keep site:mix (and therefore site and mix too):
anova(mer1, mer2) 

# Test block:
mer3 <- lmer(Yield~Diversity2 +(1+Diversity2|Site)          +(1|Mix)+(1|Site:Mix), data= biodepth) 

# block unimportant but retain block to reflect design of study:
anova(mer1, mer3) 

# Random (Diversity) intercepts (for Sites) model
mer4 <- lmer(Yield~Diversity2 +(1|Site)+(1|Block)+(1|Mix)+(1|Site:Mix), data= biodepth) 

# retain random interacepts AND slopes for diversity2|mix):
anova(mer1, mer4) 

# Extractor functions:
summary(mer1)

# latest information re. P values for MEMs
# ?pvalues

# arm package for display() and other functions:
library(arm)  
#display(mer1)

################### Box 10.4 ##################
# library(nlme)
# summary( lmer(travel~1+(1|Rail),data=Rail) )
# anova( lm(travel~1+Rail, data=Rail) )
# 16.17 + 3*615.31
# 1862.1
########

# summary() Std.Dev. are square roots of variance (components)
# sqrt(16229.04)

# BLUPs for sites, blocks and mixtures plus fixed and random slopes
coef(mer1)      

# fixed intercept and slope
fixef(mer1)

# standard errors
se.fixef(mer1)  

# random effects relative to fixed intercepts for S, B and M
ranef(mer1)     

# standard errors
se.ranef(mer1)  

# no P value
anova(mer1)     

# P values can be obtained (if you must) from: 
# languageR package, pvals.fnc() 
# lmerTest package

############
# Fig 10.5 #
############
