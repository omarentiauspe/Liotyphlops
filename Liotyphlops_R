Supplemental Material R-Code for manuscript: "More than meets no eyes: Taxonomic revision of two 
Liotyphlops Peters, 1881 (Serpentes: Anomalepididae) blindsnakes from the Atlantic Rainforest"

Omar Entiauspe-Neto et al. 2022

###Sexual dimorphism plots###

library(tidyverse)
library(palmerpenguins)
library(datasauRus)
library(gridExtra
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggplot2)
library(psych)
library(caret)
library(devtools)
library (ggbiplot)
library(ggfortify)
library(reshape2)
library(GGally)


setwd("~/R_Files")

#or setwd("~/Rdir")# 

mydata <- read.csv("Liob.csv", header = TRUE)
df <- data.frame(mydata)
str(df)


df %>%
 dplyr::select(cont_TTL, cont_TL, disc_Dorsals, disc_Ven, disc_SubC) %>%
 GGally::ggpairs(aes(color = df$Sex)) +
 scale_colour_manual(values = c("darkorange", "cyan4")) +
 scale_fill_manual(values = c("darkorange", "cyan4")) +
 theme_bw()
 
 
 #OR#
 
 df %>%
 dplyr::select(Sex, Total_length, Tail_length, Mid_dorsals, Ventrals, Subcaudals) %>%
 GGally::ggpairs(aes(color = Sex)) +
 scale_colour_manual(values = c("darkorange", "purple")) +
 scale_fill_manual(values = c("darkorange", "purple")) +
 theme_bw()
 
 #### Genetic distances ####
 
 
library(ggplot2)
library(reshape)

coxxer <- read.delim("DSS.tsv", header=TRUE)

melted_cormat <- melt(coxxer)

melted_cormat$X <- as.character(melted_cormat$X)

melted_cormat$X <- factor(melted_cormat$X, levels=unique(melted_cormat$X))

head(melted_cormat)

ggplot(data = melted_cormat, aes(x=X, y=variable, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), color = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1)) +
  theme(axis.text.y = element_text(size = 10, hjust = 1)) +
  scale_fill_gradient2(mid = "#080708", high = "#546063")
coord_fixed()

### ANOVA ###


##Hi, there! Happy to see you're interested in our supplementary material##

library(tidyverse)
library(ggpubr)
library(rstatix)

##Import data##

setwd("~/R_Files")
mydata <- read.csv("ASC.csv")
df <- data.frame(mydata)

##PLOT variables to evaluate for outliers##

ggboxplot(df, x = "Species", y = "SC")

df %>% 
  group_by(Species) %>%
  identify_outliers(SC)
  
  

##Build the linear model##
model  <- lm(TL ~ Species, data = df)

##Create a QQ plot of residuals##
ggqqplot(residuals(model))
shapiro_test(residuals(model))

##Check Normality##

df %>%
  group_by(Species) %>%
  shapiro_test(TL)

##Check homogeneity##

df %>% levene_test(TL ~ Species)

##Run ANOVA##

res.aov <- df %>% anova_test(SC ~ Species)
res.aov

##Run WILCOXON##

res <- wilcox.test(Sex, TL)
res

#### SDM ###


Liotyphlops <- read.csv("Wil.csv", header = TRUE)

Liotyphlops$lon <- Liotyphlops$lon
Liotyphlops$lat <- Liotyphlops$lat

##Clean##

capgeo <- Liotyphlops %>%
  select(Species, lat, lon) %>% 
  filter(complete.cases(.)) %>%  
  
  distinct() 

#mapa#

map(
  'worldHires',
  xlim = c(min(capgeo$lon) - 5.0, max(capgeo$lon) + 5.0),
  ylim = c(min(capgeo$lat) - 5.0, max(capgeo$lat) + 5.0),
  fill = T,
  col = "light grey"
)

box()

points(capgeo$lon,
       capgeo$lat,
       col = "orange",
       pch = 20,
       cex = 0.7)

capc <- capgeo %>% 
  dplyr::select(Species, lat, lon)

write.csv(capc, "Liotyphlops_locs.csv", row.names = FALSE)


capg <- read.csv("Liotyphlops_locs.csv")

buff <- 5   

xmin <- min(capg$lon) - 10.0
xmax <- max(capg$lon) + 10.0
ymin <- min(capg$lat) - 10.0
ymax <- max(capg$lat) + 10.0

e <- extent(xmin, xmax, ymin, ymax)

envcrop <- crop(env, e)

plot(envcrop[[1]], main = "Annual Mean Temperature")
map(
  'worldHires',
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  fill = F,
  add = T
)

points(capg$lon, capg$lat, pch = "+")


ten_div <-
  c(1, 2, 5, 6, 7, 8, 9, 10, 11)  

for (layer in ten_div) {
  envcrop[[layer]] <- envcrop[[layer]] / 10
}

layers <- c(1, 4, 11, 12, 14, 19)

writeRaster(
  stack(envcrop[[layers]]),
  paste0("bio", layers), 
  bylayer = TRUE,
  format = 'ascii',
  overwrite = T
)



##Script for PCA##

library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggplot2)
library(psych)
library(caret)
library(devtools)
library (ggbiplot)
library(ggfortify)
library(reshape2)

##Import data##

setwd("~/R_Files")
mydata <- read.csv("Males.csv", header = TRUE)
df <- data.frame(mydata)
str(df)

#This should be applied for all rows, except the "SPECIES" row#

x_num <- as.numeric(x)

[Should read as]

df$cont_SVL= as.numeric(df$cont_SVL)
df$cont_TTL= as.numeric(df$cont_TTL)
df$cont_TL= as.numeric(df$cont_TL)
df$cont_HL= as.numeric(df$cont_HL)
df$disc_DO1= as.numeric(df$disc_DO1)
df$disc_DO3= as.numeric(df$disc_DO3)
df$disc_POSTF= as.numeric(df$disc_POSTF)
df$disc_CNAS= as.numeric(df$disc_CNAS)
df$disc_Sup= as.numeric(df$disc_Sup)
df$disc_Inf= as.numeric(df$disc_Inf)
df$disc_Dorsals= as.numeric(df$disc_Dorsals)
df$disc_Ven= as.numeric(df$disc_Ven)
df$disc_SubC= as.numeric(df$disc_SubC)
df$rat_VenwSubC= as.numeric(df$rat_VenwSubC)
df$cat_CCL= as.numeric(df$cat_CCL)

#Applied to the SPECIES row#

x_fac <- as.factor(x)

[Should read as] 

df$Species = as.factor(df$Species)


# Exploratory analyses#

set.seed(111)
ind <- sample(2, nrow(df),
              replace = TRUE,
              prob = c(0.8, 0.2))
training <- df[ind==1,]
testing <- df[ind==2,]

pairs.panels(training[,-5],
             gap = 0,
             bg = c("red", "yellow", "blue")[training$Species],
             pch=21)
			 


# Species needs to be back at numeric.#

df$Species = as.numeric(as.factor(df$Species))

##Remove zero variance variables##

nearZeroVar(df)

##Plot PCA, excluding zero variance variables##

df$Species = as.factor(df$Species)

pc <- prcomp(df[-2],
              center = TRUE,
              scale. = TRUE)
			  
print(pc)

g <- ggbiplot(pc, 
               obs.scale = 1,
               var.scale = 1,
               groups = df$Species,
              ellipse = TRUE,
               circle = TRUE,
               ellipse.prob = 0.68)
> g <- g + scale_color_discrete(name = '')
> g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
> print(g)

##USING GGPLOT##

df$Species = as.numeric(as.factor(df$Species))

pca_res <- prcomp(df[-2], scale. = TRUE)

autoplot(pca_res)

df$Species = as.factor(df$Species)

##CREATE A NEW DATAFRAME#

PCi<-data.frame(pca_res$x,Species=df$Species)

##AUTOPLOT OPT 1##


autoplot(pca_res, data = df, colour = 'Species', frame = true, size = 4, scale_color_manual(values = c("#FF1BB3","#A7FF5B","#99554D")))

export = 4.5 X 4.5 


##AUTOPLOT OPT 2###

ggplot(PCi,aes(x=PC1,y=PC2,colour=Species,pch=17))+
    scale_shape_identity()+
    stat_ellipse(geom = "polygon",
                 fill = 4, alpha = 0.01)+
    geom_point(shape=17,size=4)+
    scale_color_manual(values = c("#FAA613","#F44708","#550527","#FF36AB","#A10702","#688E26"))
	
export = 3.5 x 3.5 


##MULTIPLE BOXPLOTS##


df.m <- melt(df, id.var = "Species")

p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
    geom_boxplot(aes(fill=Species))+
geom_point(aes(shape = Species), size = 2, position = position_jitterdodge(0.5), alpha=0.8)
p + facet_wrap( ~ variable, scales="free")

export = 15 x 10


####

PALETTE
SP01 [SOUSAI] = "#550527"
SP02 [BEUI] = "#FAA613"
SP03 [TERNETZII] = "#F44708"
SP04 [CAISSARA] = "#FF36AB"
SP05 [WILDERI] = "#A10702"
SP06 [TREFAUTI] = "#688E26"

