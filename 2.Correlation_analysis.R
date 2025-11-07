###############################################################################
# Create serotype-specific weights and evaluate them by correlation/regression
# on linear and logit scales
#
# This script:
#  1. Reads carriage + IPD accessory-gene matrices
#  2. Builds serotype-specific weights (from pre-vaccine years)
#  3. Joins weights to IPD presence/absence to make weighted matrices
#  4. Computes accessory-COG frequencies for 3 epochs (E1/E2/E3)
#  5. Compares:
#       - carriage vs carriage (across epochs) --> optional
#       - IPD vs IPD (before vs after weighting)
#       - IPD vs carriage (before vs after)
#  6. Repeats comparison on logit scale

# Author: Xueting Qiu
# Original: Jan 2022
# Refactored: Nov 2025
###############################################################################

# ========================== 0. Libraries =====================================

library(tidyverse)   # ggplot2, dplyr, tidyr, readr
library(readxl)
library(ggpubr)
library(gridExtra)
library(scatterD3)

# ==============Part I. Input data & generate serotype-specific rates===========

##input data for carriage pre-vaccine counts and proportions
df <- read_csv("data/CDC98_18_IPD_UScarriage-accessory-gene-with5-95percent.csv")
#this data generated from R script - "1.presence-absence-matrix.R"
ca <- df %>% filter(State %in% c("MA", "Navajo"))

#table(ca$Collection_date)

ca_pre <- subset(ca, ca$Collection_date <= 2001)
ca_pre_sero <- data.frame(table(ca_pre$WGS_serotype))

##input data for invasive data from CDC, including all the non-sequenced samples
ipd <-read_excel("data/cdc-all-metadata.xlsx")
table(ipd$Collection_Year)
ipd_pre <-subset(ipd, ipd$Collection_Year<= 1999)
ipd_pre_sero <- data.frame(table(ipd_pre$FINAL_SEROTYPE))

write.csv(ca_pre_sero,"ca_pre_sero.csv", row.names = FALSE)
write.csv(ipd_pre_sero,"ipd_pre_sero.csv", row.names = FALSE)

## get the new serotypes in the ipd after vaccines are introduced
ipd_sq <- df %>% filter(!State %in% c("MA", "Navajo"))
table(ipd_sq$WGS_serotype)
ipd_sq_sero <-data.frame(table(ipd_sq$WGS_serotype))
write.csv(ipd_sq_sero ,"ipd_sq_sero.csv", row.names = FALSE)

#########compile all three output together to compute the weight for each serotype
###using data from the pre-vaccine epoch of 1998–2001 in the United States, the weight for each serotype is the ratio of the proportion of the serotype in the carriage data (i.e., the combined MA and Southwest US data sets from 1998 to 2001) to the proportion in the IPD data from 1998 to 1999.
## details please refer to the methods section of "Serotype-specific Weight". 
######

##generate the weight distribution figure
cdc_weight_final <- read_excel("data/cdc-weight-final-022022.xlsx")

p0 <- ggplot(cdc_weight_final, aes(fill=vaccinetype, y=weightadjust, x=ipd_sero)) + 
  geom_bar(stat="identity") + scale_fill_manual("Vaccine Type", values = c("PCV7" = "darkorange", "PCV13" = "deepskyblue", "NVT" = "grey60", "NA" = "lightgrey")) +
  labs(x ="Serotype", y = "Weight value") + geom_hline(yintercept = 1, linetype="dotted", color = "red") +
  geom_hline(yintercept = 15, linetype="dotted", color = "blue") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size=12, face="bold"),axis.title.y = element_text( size=12, face="bold"), legend.title=element_text(size=12), legend.text=element_text(size=11))


p1 <- ggplot(cdc_weight_final, aes(fill = vaccinetype, 
                                   y = weightadjust, 
                                   x = ipd_sero)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("Vaccine Type",
                    values = c("PCV7"  = "darkorange",
                               "PCV13" = "deepskyblue",
                               "NVT"   = "grey60",
                               "NA"    = "lightgrey")) +
  labs(x = "Serotype", y = "Weight value") +
  geom_hline(yintercept = 1,  linetype = "dotted", color = "red") +
  geom_hline(yintercept = 15, linetype = "dotted", color = "blue") +
  scale_y_log10() +   # ← this line makes y-axis log-scaled
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

weightvalue <- grid.arrange(p0, p1, ncol = 1) 

ggsave("Fig1.weightvalue.pdf", weightvalue, height = 16, width = 19, units = "in")


# =====Part II. evaluate weight: apply weight to the presence/absence matrix ===

#start with all CDC IPD data matrix
ipd_all <- dplyr::anti_join(df,ca)

#unique(ipd_all$State)
#unique(ipd_all$WGS_serotype)
#unique(MA$WGS_serotype)


dfmeta <-  ipd_all %>% dplyr::select(1:8)
dfmeta$serotype2 <- dfmeta$WGS_serotype
##18B is NVT, 18C is pcv7, so not combined - assign 18B:18C to 18C; 7A and 7F are not combined, instead change the one called 7A:F to 7F
##change "NF" to "NT"; 
dfmeta$serotype2 <- sapply(dfmeta$serotype2, function(x){if (x =="15B" |x =="15C"|x=="15B/15C") {return("15BC")}else if (x =="24F/A/B") {return("24FAB")}else if (x =="NF") {return("NT")}else if (x =="35B:35D" | x == "35B") {return("35BD")}else if (x =="7A:7F") {return("7F")}else if (x =="18B:18C") {return("18C")}else if (x =="11B/11C"|x=="11B") {return("11BC")}else{return(x)}})


#unique(dfmeta$serotype2)
wt <- data.frame(cdc_weight_final)
#unique(wt$ipd_sero)

dfmeta_joint <- left_join(dfmeta, wt,
                          by = c("serotype2" = "ipd_sero"))

dfmeta_joint <- dfmeta_joint %>% dplyr::select(1:10, "weight","weightadjust")
df01 <- ipd_all %>% dplyr::select(accession, btuD_2:gatA_3)

df2 <- left_join(dfmeta_joint, df01,
                 by = c("accession" = "accession"))

###multiple weight to each column btuD_2:gatA_3
#library(dplyr)
df02 <- ipd_all %>% dplyr::select(btuD_2:gatA_3)
df03 <- sapply(df02, function(x) as.numeric(x)*df2$weightadjust)

df4 <- cbind(dfmeta_joint, df03)
#unique(df4$serotype2)
#unique(df2$serotype2)

# ============= part II.1 carriage correlation of COG frequencies ==============
####define the scatter plot and stats function 

corr_scatter_with_stats <- function(data, x, y,
                                    xlab = x,
                                    ylab = y,
                                    method = "pearson",
                                    digits_r = 4,
                                    digits_mse = 6,
                                    digits_rse = 6,
                                    x_pos = 0.05,
                                    y_pos = 1.00) {
  # grab the vectors
  xvec <- data[[x]]
  yvec <- data[[y]]
  
  # correlation test
  cor_test <- cor.test(xvec, yvec, method = method)
  r   <- unname(cor_test$estimate)
  p   <- cor_test$p.value
  
  # linear model (for MSE/RSE)
  lm_model <- lm(yvec ~ xvec)
  mse <- mean(residuals(lm_model)^2)
  rse <- summary(lm_model)$sigma
  
  # build label
  label_text <- paste0(
    "R = ", round(r, 3),
    "\np ", ifelse(p < 2.2e-16,
                   "< 2.2e-16",
                   format(p, digits = 10, scientific = TRUE)),
    "\nMSE = ", round(mse, 3),
    "\nRSE = ", round(rse, 3)
  )
  
  # compute positions in data space
  x_loc <- x_pos * (max(xvec, na.rm = TRUE) - min(xvec, na.rm = TRUE)) + min(xvec, na.rm = TRUE)
  y_loc <- y_pos * (max(yvec, na.rm = TRUE) - min(yvec, na.rm = TRUE)) + min(yvec, na.rm = TRUE)
  
  # plot
  p_out <- ggscatter(data, x = x, y = y,
                     add = "reg.line", conf.int = TRUE,
                     xlab = xlab, ylab = ylab) +
    annotate("text",
             x = x_loc, y = y_loc,
             label = label_text,
             hjust = 0, vjust = 1,
             size = 4,
             color = "black")
  
  return(p_out)
}

#calculate carriage
COGs <- ncol(ca)
table(ca$Epoch1)
# E1  E2  E3 
#407 928 916

table(ca$Epoch2)
#E1A E1B E1C E2A E2B E2C E2D E3A E3B E3C 
#105 131 171 203  79 577  69 119 478 319 
epoch1 <- as.matrix(ca[ca$Epoch1 == "E1",][,9:COGs])
epoch2 <- as.matrix(ca[ca$Epoch1 == "E2",][,9:COGs])
epoch3 <- as.matrix(ca[ca$Epoch1 == "E3",][,9:COGs])# no epoch3 in SPARC1 dataset

epoch1 <- as.data.frame(apply(epoch1, 2, as.numeric))
epoch2 <- as.data.frame(apply(epoch2, 2, as.numeric))
epoch3 <- as.data.frame(apply(epoch3, 2, as.numeric))

freqs_ca <- matrix(nrow = (COGs-8), ncol = 3)  #setting up matrix for all frequency results
freqs_ca[,1]<-as.matrix(colSums(epoch1)/nrow(epoch1))
freqs_ca[,2]<-as.matrix(colSums(epoch2)/nrow(epoch2))
freqs_ca[,3]<-as.matrix(colSums(epoch3)/nrow(epoch3))

freqs_ca <- data.frame(freqs_ca)
head(freqs_ca)


#exploration, optional

p2 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X1",
  y = "X2",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)



p3 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X1",
  y = "X3",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


p4 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X2",
  y = "X3",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


f2optional <- ggarrange(p2, p3, p4, ncol=3, labels = c("A","B","C"))
f2optional
#ggsave("OptionalFigure.pdf", f2optional, width = 15,height = 5)


# ============ Part II.2 Comparing invasive at different epochs ================

#frequency before weighing

COGs <- ncol(df2)

epoch1 <- as.matrix(df2[df2$Epoch1 == "E1",][,13:COGs])
epoch2 <- as.matrix(df2[df2$Epoch1 == "E2",][,13:COGs])
epoch3 <- as.matrix(df2[df2$Epoch1 == "E3",][,13:COGs])

epoch1 <- as.data.frame(apply(epoch1, 2, as.numeric))
epoch2 <- as.data.frame(apply(epoch2, 2, as.numeric))
epoch3 <- as.data.frame(apply(epoch3, 2, as.numeric))

freqs_ipd_be <- matrix(nrow = (COGs-12), ncol = 3)  #setting up matrix for all frequency results
freqs_ipd_be[,1]<-as.matrix(colSums(epoch1)/nrow(epoch1))
freqs_ipd_be[,2]<-as.matrix(colSums(epoch2)/nrow(epoch2))
freqs_ipd_be[,3]<-as.matrix(colSums(epoch3)/nrow(epoch3))

freqs_ipd_be <- data.frame(freqs_ipd_be)
head(freqs_ipd_be)

#####make frequency comparison plot
p5 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X1",
  y = "X2",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)

p6 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X1",
  y = "X3",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

p7 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X2",
  y = "X3",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

f3a <- ggarrange(p5, p6, p7, ncol=3, labels = c("A","B","C"))

#ggsave("figure2b.pdf", f2b, width = 15,height = 5)


###Calculate frequencies for 3 epochs - after weighing

COGs <- ncol(df4)

epoch1 <- as.matrix(df4[df4$Epoch1 == "E1",][,13:COGs])
epoch2 <- as.matrix(df4[df4$Epoch1 == "E2",][,13:COGs])
epoch3 <- as.matrix(df4[df4$Epoch1 == "E3",][,13:COGs])

wtepoch1 <- df4[df4$Epoch1 == "E1",][,"weightadjust"]
wtepoch2 <- df4[df4$Epoch1 == "E2",][,"weightadjust"]
wtepoch3 <- df4[df4$Epoch1 == "E3",][,"weightadjust"]


freqs_ipd_af <- matrix(nrow = (COGs-12), ncol = 3)  #setting up matrix for all frequency results
freqs_ipd_af[,1]<-as.matrix(colSums(epoch1)/sum(wtepoch1))
freqs_ipd_af[,2]<-as.matrix(colSums(epoch2)/sum(wtepoch2))
freqs_ipd_af[,3]<-as.matrix(colSums(epoch3)/sum(wtepoch3))

freqs_ipd_af <- data.frame(freqs_ipd_af)
head(freqs_ipd_af)

p8 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X1",
  y = "X2",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)

p9 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X1",
  y = "X3",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


p10 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X2",
  y = "X3",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


f3b <- ggarrange(p8, p9, p10, ncol=3, labels = c("D","E","F"))
f3 <- ggarrange(f3a, f3b, ncol =1)

ggsave("Fig3.pdf", f3, width = 20,height = 12)


# ======== Part II.3. Comparing carriage and invasive before and after =========
###combine before and after weighing, and carriage

freqs_ipd_af <- dplyr::rename(freqs_ipd_af, c(af_epoch1 = X1, af_epoch2 = X2, af_epoch3 = X3))
freqs_ipd_be <-  dplyr::rename(freqs_ipd_be, c(be_epoch1 = X1, be_epoch2 = X2, be_epoch3 = X3))
freqs_ca <- dplyr::rename(freqs_ca, c(ca_epoch1 = X1, ca_epoch2 = X2, ca_epoch3=X3))

freqs_comb <- cbind.data.frame(freqs_ipd_af, freqs_ipd_be, freqs_ca)
freqs_comb <- freqs_comb[,c("af_epoch1", "af_epoch2", "af_epoch3", "be_epoch1", "be_epoch2", "be_epoch3", "ca_epoch1", "ca_epoch2", "ca_epoch3")]
head(freqs_comb)

#prevaccine, before weighting
p20 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "be_epoch1",
  y = "ca_epoch1",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#prevaccine, after weighting
p21 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "af_epoch1",
  y = "ca_epoch1",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV7, before weighting
p22 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "be_epoch2",
  y = "ca_epoch2",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV7, after weighting
p23 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "af_epoch2",
  y = "ca_epoch2",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV13, before weighting
p24 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "be_epoch3",
  y = "ca_epoch3",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV13, after weighting
p25 <- corr_scatter_with_stats(
  data = freqs_comb,
  x = "af_epoch3",
  y = "ca_epoch3",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

f2a <- ggarrange(p20, p22, p24, ncol=3, labels = c("A","B","C"))
f2b <- ggarrange(p21, p23, p25, ncol=3, labels = c("D","E","F"))

f2 <- ggarrange(f2a, f2b, ncol=1)

ggsave("Fig2.pdf", f2, width = 20,height = 12)


# ======================= Part II.4 supplemental - logit transform =============

###logit transform accessory gene frequencies
#because it makes the relationship between pre- and post-selection frequencies more linear and variance-stable, allowing fair, interpretable comparisons — especially when genes range from rare to nearly fixed.
#for all carriage; for all invasive before and after

##carriage 
freqs_ca$X1logit <- log(freqs_ca$ca_epoch1/(1-freqs_ca$ca_epoch1))
freqs_ca$X2logit <- log(freqs_ca$ca_epoch2/(1-freqs_ca$ca_epoch2))
freqs_ca$X3logit <- log(freqs_ca$ca_epoch3/(1-freqs_ca$ca_epoch3))
freqs_ca <- freqs_ca[!is.infinite(rowSums(freqs_ca)),] #remove inf values

##plots are optional for carriage 
p11 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X1logit",
  y = "X2logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)

p12 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X1logit",
  y = "X3logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

p13 <- corr_scatter_with_stats(
  data = freqs_ca,
  x = "X2logit",
  y = "X3logit",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


#fs1 <- ggarrange(p11, p12, p13, ncol=3, labels = c("A","B","C"))
#fs1
#ggsave("figureS1.pdf", fs1, width = 15,height = 5)

###IPD before weight

freqs_ipd_be$X1logit <- log(freqs_ipd_be$be_epoch1/(1-freqs_ipd_be$be_epoch1))
freqs_ipd_be$X2logit <- log(freqs_ipd_be$be_epoch2/(1-freqs_ipd_be$be_epoch2))
freqs_ipd_be$X3logit <- log(freqs_ipd_be$be_epoch3/(1-freqs_ipd_be$be_epoch3))
freqs_ipd_be <- freqs_ipd_be[!is.infinite(rowSums(freqs_ipd_be)),] #remove inf values

p17 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X1logit",
  y = "X2logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)

p18 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X1logit",
  y = "X3logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

p19 <- corr_scatter_with_stats(
  data = freqs_ipd_be,
  x = "X2logit",
  y = "X3logit",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

fs2a <- ggarrange(p17, p18, p19, ncol=3, labels = c("A","B","C"))


###IPD after weight 

freqs_ipd_af$X1logit <- log(freqs_ipd_af$af_epoch1/(1-freqs_ipd_af$af_epoch1))
freqs_ipd_af$X2logit <- log(freqs_ipd_af$af_epoch2/(1-freqs_ipd_af$af_epoch2))
freqs_ipd_af$X3logit <- log(freqs_ipd_af$af_epoch3/(1-freqs_ipd_af$af_epoch3))
freqs_ipd_af <- freqs_ipd_af[!is.infinite(rowSums(freqs_ipd_af)),]

p14 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X1logit",
  y = "X2logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV7 COG frequency"
)

p15 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X1logit",
  y = "X3logit",
  xlab = "Pre-vaccine COG frequency",
  ylab = "Post-PCV13 COG frequency"
)

p16 <- corr_scatter_with_stats(
  data = freqs_ipd_af,
  x = "X2logit",
  y = "X3logit",
  xlab = "Post-PCV7 COG frequency",
  ylab = "Post-PCV13 COG frequency"
)


fs2b <- ggarrange(p14, p15, p16, ncol=3, labels = c("D","E","F"))

fs2 <- ggarrange(fs2a, fs2b, ncol=1)

ggsave("FigS2.pdf", fs2, width = 20,height = 12)



######logit scale for comparing carriage and ipd
logit <- function(p) log(p / (1 - p))

freqs_comb_logit <- freqs_comb %>%
  mutate(across(everything(),
                ~ logit(.x),
                .names = "{.col}_logit"))
freqs_comb_logit <- freqs_comb_logit[!is.infinite(rowSums(freqs_comb_logit)),]
head(freqs_comb_logit)


#prevaccine, before weighting
p26 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "be_epoch1_logit",
  y = "ca_epoch1_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#prevaccine, after weighting
p27 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "af_epoch1_logit",
  y = "ca_epoch1_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV7, before weighting
p28 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "be_epoch2_logit",
  y = "ca_epoch2_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV7, after weighting
p29 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "af_epoch2_logit",
  y = "ca_epoch2_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV13, before weighting
p30 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "be_epoch3_logit",
  y = "ca_epoch3_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

#post-PCV13, after weighting
p31 <- corr_scatter_with_stats(
  data = freqs_comb_logit,
  x = "af_epoch3_logit",
  y = "ca_epoch3_logit",
  xlab = "Invasive COG frequency",
  ylab = "Carriage COG frequency"
)

fs1a <- ggarrange(p26, p28, p30, ncol=3, labels = c("A","B","C"))
fs1b <- ggarrange(p27, p29, p31, ncol=3, labels = c("D","E","F"))

fs1 <- ggarrange(fs1a, fs1b, ncol = 1, nrow = 2) 

ggsave("FigS1.pdf", fs1, width = 20,height = 12)




