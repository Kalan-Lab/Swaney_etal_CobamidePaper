
# Script 6 ----------------------------------------------------------------
# Figure 4

# Cobamide Skin Microbiome Manuscript
# C. amycolatum LK19 cell extract in E. coli cobamide assay
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

#read in assay results
data <- read_excel("6_EcoliAssay/data/LK19CellExtracts.xlsx")


# Plotting ----------------------------------------------------------------

# order the factors for plotting - concentrations/dilutions
data$Concentration <- factor(data$Concentration, levels= c("2.5 ng/mL", "2 ng/mL","1.5 ng/mL","1 ng/mL", "0.5 ng/mL",
                                                           "0.25 ng/mL","0.1 ng/mL","0 ng/mL",
                                                           "1:10000","1:25000", "1:37500","1:50000"))

# order the factors for plotting - sample order
data$Sample <- factor(data$Sample, levels = c("CnCbl","LK19"))

# exclude 2.5 and 2 ng/mL standards because they are outside of linear range
data <- data %>% filter(Concentration != "2.5 ng/mL" & Concentration != "2 ng/mL")

# boxplot of log10 OD600 vs sample concentration/dilution
p1 <- ggplot(data, aes(x=Concentration, y=log10(OD600), fill=Sample)) + geom_boxplot() + 
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1)) + 
  ylab("Log10 OD600") + geom_point()

# standard curve for cyanocobalamin standards
p2 <- ggplot(data, aes(x=Num_Concentration,y=log10(OD600))) + geom_point() +
  geom_smooth(method = "lm") + xlab("Cyanocobalamin concentration (ng/mL)")

# plot - Figure 4
plot_grid(p2,p1, align="h",axis = "b")

# get equation for line
lm(log10(OD600) ~ Num_Concentration, data)
# y = 0.6118x - 1.2585


# Calculate cobamide concentration in extracts -------------------------------

# take the mean and get the standard deviation of each CE dilution or CnCbl concentration replicate
data <- data %>% group_by(Concentration,Sample,Label) %>% summarise(logmean=mean(log10(OD600)))

# add in wet cell weight information for each cell extract
wet_cell_weights <- data.frame(wcw = c(10.72, 9.91, 11.53, 10.67, 9.42, 9.51), Label = c(3, 4, 5, 6, 7, 8))
wet_cell_weights$Label <- as.character(wet_cell_weights$Label)
data <- left_join(data,wet_cell_weights)

# add in OD600 information for each cell extract
OD <- data.frame(OD600 = c(1.21, 1.57, 1.35, 1.12, 1.77, 0.92), Label = c(3, 4, 5, 6, 7, 8))
OD$Label <- as.character(OD$Label)
data <- left_join(data,OD)

# calculate cobamide concentration in the extract dilutions
# y = -1.2585x + 0.6118
CE_10000 <- data %>% filter(Concentration == "1:10000")
CE_25000 <- data %>% filter(Concentration == "1:25000")
CE_37500 <- data %>% filter(Concentration == "1:37500")
CE_50000 <- data %>% filter(Concentration == "1:50000")

# 1:10000 dilution
CE_10000$cobamide_concentration_ng_mL <- (CE_10000$logmean + 1.2585)/0.6118
CE_10000$DilutionFactor <- 10000
CE_10000$total_cobamides <- CE_10000$cobamide_concentration_ng_mL * CE_10000$DilutionFactor * 1.1 /1000
CE_10000$ug_cobamides_per_wcw <- CE_10000$total_cobamides/CE_10000$wcw
CE_10000$ug_cobamides_per_wcw <- CE_10000$cobamide_concentration_ng_mL * CE_10000$DilutionFactor * 1.1 /1000 /CE_10000$wcw
CE_10000$ug_cobamides_per_L <- CE_10000$cobamide_concentration_ng_mL * CE_10000$DilutionFactor * 1.1 /1000 / 1 # 1 L culture
CE_10000$total_cells <- CE_10000$OD600 * 8E8 * 1000 # OD600 x 8E8 cells/mL (approximate cell number in OD600 of 1) x 1000 mL
CE_10000$total_cell_volume_cm3 <- CE_10000$total_cells / (1/1E-12) # 1 cell = 1 um3 and 1 um3 = 1E-12 cm3
CE_10000$mol_cobamides <- CE_10000$total_cobamides / 1E6 / 1355 # 1 ug = 1E-6 g and CnCbl mol weight = 1355 g/mol
CE_10000$intracellular_cobamide_concentration <- CE_10000$mol_cobamides / (CE_10000$total_cell_volume_cm3 / 1000) * 1E6 # uM

# 1:25000 dilution
CE_25000$cobamide_concentration_ng_mL <- (CE_25000$logmean + 1.2585)/0.6118
CE_25000$DilutionFactor <- 25000
CE_25000$total_cobamides <- CE_25000$cobamide_concentration_ng_mL * CE_25000$DilutionFactor * 1.1 /1000
CE_25000$ug_cobamides_per_wcw <- CE_25000$total_cobamides/CE_25000$wcw
CE_25000$ug_cobamides_per_wcw <- CE_25000$cobamide_concentration_ng_mL * CE_25000$DilutionFactor * 1.1 /1000 /CE_25000$wcw
CE_25000$ug_cobamides_per_L <- CE_25000$cobamide_concentration_ng_mL * CE_25000$DilutionFactor * 1.1 /1000 / 1 # 1 L culture
CE_25000$total_cells <- CE_25000$OD600 * 8E8 * 1000 # OD600 x 8E8 cells/mL (approximate cell number in OD600 of 1) x 1000 mL
CE_25000$total_cell_volume_cm3 <- CE_25000$total_cells / (1/1E-12) # 1 cell = 1 um3 and 1 um3 = 1E-12 cm3
CE_25000$mol_cobamides <- CE_25000$total_cobamides / 1E6 / 1355 # 1 ug = 1E-6 g and CnCbl mol weight = 1355 g/mol
CE_25000$intracellular_cobamide_concentration <- CE_25000$mol_cobamides / (CE_25000$total_cell_volume_cm3 / 1000) * 1E6 # uM

# 1:37500 dilution
CE_37500$cobamide_concentration_ng_mL <- (CE_37500$logmean + 1.2585)/0.6118
CE_37500$DilutionFactor <- 37500
CE_37500$total_cobamides <- CE_37500$cobamide_concentration_ng_mL * CE_37500$DilutionFactor * 1.1 /1000
CE_37500$ug_cobamides_per_wcw <- CE_37500$total_cobamides/CE_37500$wcw
CE_37500$total_cobamides <- CE_37500$cobamide_concentration_ng_mL * CE_37500$DilutionFactor * 1.1 /1000
CE_37500$ug_cobamides_per_wcw <- CE_37500$total_cobamides/CE_37500$wcw
CE_37500$ug_cobamides_per_L <- CE_37500$cobamide_concentration_ng_mL * CE_37500$DilutionFactor * 1.1 /1000 / 1 # 1 L culture
CE_37500$total_cells <- CE_37500$OD600 * 8E8 * 1000 # OD600 x 8E8 cells/mL (approximate cell number in OD600 of 1) x 1000 mL
CE_37500$total_cell_volume_cm3 <- CE_37500$total_cells / (1/1E-12) # 1 cell = 1 um3 and 1 um3 = 1E-12 cm3
CE_37500$mol_cobamides <- CE_37500$total_cobamides / 1E6 / 1355 # 1 ug = 1E-6 g and CnCbl mol weight = 1355 g/mol
CE_37500$intracellular_cobamide_concentration <- CE_37500$mol_cobamides / (CE_37500$total_cell_volume_cm3 / 1000) * 1E6 # uM

# 1:50000 dilution
CE_50000$cobamide_concentration_ng_mL <- (CE_50000$logmean + 1.2585)/0.6118
CE_50000$DilutionFactor <- 50000
CE_50000$total_cobamides <- CE_50000$cobamide_concentration_ng_mL * CE_50000$DilutionFactor * 1.1 /1000
CE_50000$ug_cobamides_per_wcw <- CE_50000$total_cobamides/CE_50000$wcw
CE_50000$ug_cobamides_per_wcw <- CE_50000$cobamide_concentration_ng_mL * CE_50000$DilutionFactor * 1.1 /1000 /CE_50000$wcw
CE_50000$ug_cobamides_per_L <- CE_50000$cobamide_concentration_ng_mL * CE_50000$DilutionFactor * 1.1 /1000 / 1 # 1 L culture
CE_50000$total_cells <- CE_50000$OD600 * 8E8 * 1000 # OD600 x 8E8 cells/mL (approximate cell number in OD600 of 1) x 1000 mL
CE_50000$total_cell_volume_cm3 <- CE_50000$total_cells / (1/1E-12) # 1 cell = 1 um3 and 1 um3 = 1E-12 cm3
CE_50000$mol_cobamides <- CE_50000$total_cobamides / 1E6 / 1355 # 1 ug = 1E-6 g and CnCbl mol weight = 1355 g/mol
CE_50000$intracellular_cobamide_concentration <- CE_50000$mol_cobamides / (CE_50000$total_cell_volume_cm3 / 1000) * 1E6 # uM

CE_calculations <- rbind(CE_10000, CE_25000)
CE_calculations <- rbind(CE_calculations, CE_37500)
CE_calculations <- rbind(CE_calculations, CE_50000)

# ug cobamides per wet cell weight 
calculations_final <- CE_calculations %>% group_by(Label) %>% summarise(mean=mean(ug_cobamides_per_wcw))
calculations_final %>% summarise(mean=mean(calculations_final$mean), sd = sd(calculations_final$mean))

# uM cobamides intracellular
intracellular_final <- CE_calculations %>% group_by(Label) %>% summarise(mean=mean(intracellular_cobamide_concentration))
intracellular_final %>% summarise(mean=mean(intracellular_final$mean), sd = sd(intracellular_final$mean))
