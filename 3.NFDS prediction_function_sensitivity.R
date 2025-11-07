###############################################################################
# Script: NFDS prediction for weighted population and accessory COGs 
#
# Description:
#   This script applies serotype-specific weights to CDC IPD isolates,
#   aggregates them by GPSC-defined strain clusters, and then uses an
#   NFDS (negative frequency–dependent selection) prediction framework
#   to estimate post-PCV13 strain proportions from pre-vaccine accessory-gene
#   profiles.
#
#   Workflow:
#     1. Read in the combined IPD + carriage accessory-gene matrix.
#     2. Separate out IPD isolates and attach cleaned serotype, year/epoch,
#        and CDC-derived serotype weights.
#     3. Join GPSC (poppunk) cluster assignments to each isolate.
#     4. Summarize strain (GPSC) proportions by epoch:
#           - unweighted 
#           - weighted   
#        and make barplots (plot4a/plot4b --> Figure 4).
#     5. Build strain-level accessory-gene frequency matrices (per GPSC),
#        derive the expected pre-vaccine COG distribution, and run QP to get
#        predicted post-vaccine strain proportions.
#     6. Compare predicted vs. observed post-PCV13 strain proportions and
#        plot them (Fig 5. regression/1:1 plots) with goodness-of-fit
#        metrics (SSE, RMSE, adj. R²).
#
#   Sensitivity analyses included:
#     - Varying the minimum GPSC size used in the prediction:
#         * main analysis: GPSC with > 5 isolates
#         * sensitivity:   GPSC with > 3 isolates
#         * sensitivity:   GPSC with > 10 isolates
#     - Reclassifying serotype 3 as “non-vaccine type” (NVT) in the
#       vaccinetype variable, to reflect epidemiologic evidence that
#       PCV13 has limited impact on serotype 3.
#
# Author:   Xueting Qiu
# Original:  July 2022
# Refactored: Nov 2025
###############################################################################

###############################################################################
# 0. libraries
###############################################################################
library(tidyverse)
library(readr)
library(readxl)
library(ggrepel)
library(ggpubr)
library(quadprog)
library(Metrics)
library(EnvStats)

###############################################################################
# 1. READ DATA ONCE
###############################################################################

# accessory gene matrix (IPD + carriage)
acc_all  <- read_csv("data/CDC98_18_IPD_UScarriage-accessory-gene-with5-95percent.csv")
#this data generated from R script - "1.presence-absence-matrix.R"

# serotype-specific weights
wt       <- read_excel("data/cdc-weight-final-022022.xlsx")

# GPSC assignments
clusters <- read_csv("data/cdcIPD-GPSC-assigned_clusters2.csv")

# keep only IPD
ipd_all <- acc_all %>%
  filter(!State %in% c("MA", "Navajo")) %>%
  filter(Collection_date != "1998",
         Collection_date != "1999")

# metadata + your serotype recode (do NOT convert serotype 3 here)
ipd_meta <- ipd_all %>%
  dplyr::select(1:8) %>%
  mutate(
    serotype2 = case_when(
      WGS_serotype %in% c("15B","15C","15B/15C") ~ "15BC",
      WGS_serotype == "24F/A/B"                 ~ "24FAB",
      WGS_serotype == "NF"                      ~ "NT",
      WGS_serotype %in% c("35B:35D","35B")      ~ "35BD",
      WGS_serotype == "7A:7F"                   ~ "7F",
      WGS_serotype == "18B:18C"                 ~ "18C",
      WGS_serotype %in% c("11B/11C","11B")      ~ "11BC",
      is.na(WGS_serotype)                       ~ "NT",
      TRUE                                      ~ WGS_serotype    # <-- serotype 3 stays "3"
    ),
    Epoch2 = case_when(
      Collection_date == "2009" ~ "E1",
      Collection_date %in% c("2013","2015","2016") ~ "E2",
      TRUE ~ "E3"
    )
  )

# add weight columns
ipd_meta_wt <- ipd_meta %>%
  left_join(wt, by = c("serotype2" = "ipd_sero")) %>%
  dplyr::select(1:11, weight, weightadjust)

# unweighted PA
ipd_pa_unw <- ipd_all %>%
  dplyr::select(accession, btuD_2:gatA_3)

# combine → master unweighted
df2_master <- ipd_meta_wt %>%
  left_join(ipd_pa_unw, by = "accession")

# weighted PA
pa_only     <- ipd_pa_unw %>% dplyr::select(-accession)
pa_weighted <- sapply(pa_only, function(x) as.numeric(x) * df2_master$weightadjust)

# combine → master weighted
df4_master <- cbind(ipd_meta_wt, pa_weighted)

# join GPSC once
df2_master <- clusters %>%
  inner_join(df2_master, by = c("Taxon" = "accession"))

df4_master <- clusters %>%
  inner_join(df4_master, by = c("Taxon" = "accession"))

# unify a couple of GPSC names once
df2_master$Cluster[df2_master$Cluster == "3_6_409_1107"] <- "3"
df2_master$Cluster[df2_master$Cluster == "388_1033"]     <- "388"
df4_master$Cluster[df4_master$Cluster == "3_6_409_1107"] <- "3"
df4_master$Cluster[df4_master$Cluster == "388_1033"]     <- "388"

# helper for plotting
collapse_vt <- function(x) {
  if (is.na(x)) "NVT" else if (x %in% c("PCV7","PCV13")) "VT" else "NVT"
}

###############################################################################
# 2. BUILD PLOT4A / PLOT4B ONCE (GPSC > 5)
###############################################################################

## unweighted
gpsc_counts_unw <- table(df2_master$Cluster)
keep_5_unw      <- names(gpsc_counts_unw)[gpsc_counts_unw > 5]

df_plot_unw <- df2_master %>%
  filter(Cluster %in% keep_5_unw) %>%
  mutate(vaccinetype = sapply(vaccinetype, collapse_vt))

dfF_unw <- df_plot_unw %>%
  dplyr::select(Cluster, Epoch2) %>%
  group_by(Epoch2) %>%
  count(Cluster) %>%
  mutate(freq = prop.table(n)) %>%
  ungroup() %>%
  dplyr::select(Epoch2, Cluster, freq) %>%
  tidyr::pivot_wider(names_from = Epoch2, values_from = freq, values_fill = 0) %>%
  arrange(Cluster)

vaccineT_plot <- df_plot_unw %>%
  distinct(Cluster, vaccinetype) %>%
  subset(vaccinetype == "VT") %>%
  dplyr::rename(W = vaccinetype) %>%
  full_join(
    df_plot_unw %>% distinct(Cluster, vaccinetype) %>% subset(vaccinetype == "NVT"),
    by = "Cluster"
  ) %>%
  tidyr::unite("vaccine", W:vaccinetype, na.rm = TRUE)

dfF_unw <- dfF_unw %>% left_join(vaccineT_plot, by = "Cluster")

datPlotA <- dfF_unw %>%
  dplyr::select(SC = Cluster, E1, E3, vaccine) %>%
  pivot_longer(c(E1, E3), names_to = "Epoch", values_to = "Prevalence") %>%
  arrange(SC, Epoch) %>%
  mutate(SC = factor(SC, levels = unique(SC)))

plot4a <- ggplot(datPlotA,
                 aes(x = SC, y = Prevalence, alpha = Epoch, fill = vaccine)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_alpha_manual(values = c("E1" = 1, "E3" = 0.4),
                     labels = c("E1" = "Pre-vaccine", "E3" = "Post-vaccine"),
                     name   = "Epoch") +
  scale_fill_manual(values = c("NVT"     = "#143c77",
                               "VT"      = "darkred",
                               "VT_NVT"  = "darkorchid4"),
                    labels = c("NVT"    = "Non-vaccine type",
                               "VT"     = "Vaccine type",
                               "VT_NVT" = "Mixed"),
                    name   = "Composition") +
  xlab("Strain (SC)") + ylab("Strain proportion") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.justification = c(1, 1),
        legend.box = "horizontal",
        legend.position = c(1, 1),
        legend.spacing.y = unit(0.2, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = gray(0.96), color = NA))+
  guides(alpha = guide_legend(order = 1),
         fill  = guide_legend(order = 2))

## weighted
gpsc_counts_w <- table(df4_master$Cluster)
keep_5_w      <- names(gpsc_counts_w)[gpsc_counts_w > 5]

df_plot_w <- df4_master %>%
  filter(Cluster %in% keep_5_w) %>%
  mutate(vaccinetype = sapply(vaccinetype, collapse_vt))

dfF_w <- df_plot_w %>%
  dplyr::select(Cluster, Epoch2, weightadjust) %>%
  group_by(Epoch2, Cluster) %>%
  summarise(cnt = sum(weightadjust), .groups = "drop_last") %>%
  mutate(freq = cnt / sum(cnt)) %>%
  ungroup() %>%
  dplyr::select(Epoch2, Cluster, freq) %>%
  tidyr::pivot_wider(names_from = Epoch2, values_from = freq, values_fill = 0) %>%
  arrange(Cluster) %>%
  left_join(vaccineT_plot, by = "Cluster")

datPlotB <- dfF_w %>%
  dplyr::select(SC = Cluster, E1, E3, vaccine) %>%
  pivot_longer(c(E1, E3), names_to = "Epoch", values_to = "Prevalence") %>%
  arrange(SC, Epoch) %>%
  mutate(SC = factor(SC, levels = unique(SC)))

plot4b <- ggplot(datPlotB,
                 aes(x = SC, y = Prevalence, alpha = Epoch, fill = vaccine)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_alpha_manual(values = c("E1" = 1, "E3" = 0.4),
                     labels = c("E1" = "Pre-vaccine", "E3" = "Post-vaccine"),
                     name   = "Epoch") +
  scale_fill_manual(values = c("NVT"     = "#143c77",
                               "VT"      = "darkred",
                               "VT_NVT"  = "darkorchid4"),
                    labels = c("NVT"    = "Non-vaccine type",
                               "VT"     = "Vaccine type",
                               "VT_NVT" = "Mixed"),
                    name   = "Composition") +
  xlab("Strain (SC)") + ylab("Strain proportion") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.justification = c(1, 1),
        legend.box = "horizontal",
        legend.position = c(1, 1),
        legend.spacing.y = unit(0.2, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = gray(0.96), color = NA))+
  guides(alpha = guide_legend(order = 1),
         fill  = guide_legend(order = 2))+
  theme(legend.position = "none")

fig4 <- ggarrange(plot4a, plot4b, ncol = 1) 
ggsave("Fig4.pdf", fig4, width = 20,height = 12)


###############################################################################
# 3. QP helper
###############################################################################
QP <- function(X, Y) {
  rinv <- solve(chol(t(X) %*% X))
  C    <- cbind(rep(1, ncol(X)), diag(ncol(X)))
  b    <- c(1, rep(0, ncol(X)))
  d    <- t(Y) %*% X
  sol  <- solve.QP(Dmat = rinv,
                   factorized = TRUE,
                   dvec = d,
                   Amat = C,
                   bvec = b,
                   meq  = 1)$solution
  round(sol, 5)
}

###############################################################################
# 4. PREDICTION / SENSITIVITY FUNCTION
#    This keeps your original prediction code; only the setup is parameterized.
###############################################################################
run_prediction_sensitivity <- function(df_unw_master,
                                       df_wt_master,
                                       gpsc_cutoff = 5,
                                       sero3_as_nvt = FALSE,
                                       suffix = "sens") {
  
  # ---- filter by GPSC size (unweighted & weighted) ----
  tab_unw <- table(df_unw_master$Cluster)
  keep    <- names(tab_unw)[tab_unw > gpsc_cutoff]
  
  df  <- df_unw_master %>% filter(Cluster %in% keep)
  df4 <- df_wt_master %>% filter(Cluster %in% keep)
  
  # ---- collapse vaccinetype the same way you do in your code ----
  df$vaccinetype <- sapply(df$vaccinetype, function(x){
    if (is.na(x))           "NVT"
    else if (x %in% c("PCV7","PCV13")) "VT"
    else                    "NVT"
  })
  df4$vaccinetype <- sapply(df4$vaccinetype, function(x){
    if (is.na(x))           "NVT"
    else if (x %in% c("PCV7","PCV13")) "VT"
    else                    "NVT"
  })
  
  # ---- sensitivity: only change vaccinetype for serotype 3 ----
  if (sero3_as_nvt) {
    message("Applying serotype 3 sensitivity: converting serotype 3 → NVT and removing from Cluster 1")
    
    # convert serotype 3 to NVT before collapsing
    df$vaccinetype[df$serotype2 == "3"]  <- "NVT"
    df4$vaccinetype[df4$serotype2 == "3"] <- "NVT"
    
    # remove serotype 3 from Cluster 1 (keep rest of Cluster 1)
    df <- df %>%
      filter(!(Cluster == "1" & serotype2 == "3"))
    df4 <- df4 %>%
      filter(!(Cluster == "1" & serotype2 == "3"))
  }
  
  #############
  # PREDICTION CODE: UNWEIGHTED
  #############
  
  # cluster × vaccinetype × epoch
  dfFVT <- df %>% dplyr::select(Cluster,vaccinetype, Epoch2) %>% 
    group_by(Epoch2) %>%  count(Cluster,vaccinetype) %>% 
    mutate(freq = round(prop.table(n), digits = 6)) %>% ungroup() %>%
    dplyr::select(Epoch2, Cluster,vaccinetype, freq) %>% 
    spread(Epoch2, freq, fill = 0) %>% arrange(Cluster,vaccinetype)
  
  # accessory gene frequencies in each strain
  SC_freq_df <- df %>% dplyr::select(Cluster, vaccinetype, Epoch2, 
                                     btuD_2:gatA_3) %>%
    arrange(Cluster) %>% group_by(Cluster,vaccinetype,Epoch2) %>%
    mutate(SC_n = n()) %>% ungroup() %>% 
    group_by(Cluster,vaccinetype,Epoch2,SC_n) %>%
    summarise_at(vars(btuD_2:gatA_3),mean) %>% 
    ungroup()
  
  #### Present at E1 ####
  SCE1 <- dfFVT %>% subset(E1 > 0) %>% 
    dplyr::select(Cluster,vaccinetype) %>% 
    mutate(Epoch2 = "E1")
  
  #### NVTs Present at E1 ####
  SCE2 <- SCE1 %>% subset(vaccinetype == "NVT")
  
  ## 1)SC not present in prevaccine n = ...
  dfImputed <- dfFVT %>% subset(vaccinetype == "NVT" & (E1 == 0)) %>%
    dplyr::select(Cluster,vaccinetype) %>% mutate(Epoch2 = "E1", n=1) %>%
    dplyr::select(Epoch2,Cluster,vaccinetype,n)
  
  dfFVTImputed <- df %>% dplyr::select(Cluster,vaccinetype, Epoch2) %>% 
    group_by(Epoch2) %>%  count(Cluster,vaccinetype) %>% 
    ungroup() %>% bind_rows(dfImputed) %>% group_by(Epoch2) %>%
    mutate(freq = round(prop.table(n), digits = 6)) %>% ungroup() %>%
    dplyr::select(Epoch2, Cluster,vaccinetype, freq) %>% 
    spread(Epoch2, freq, fill = 0) %>% arrange(Cluster,vaccinetype) 
  
  dfFNVTImputed <- dfFVTImputed %>% 
    subset(vaccinetype == "NVT") %>%
    mutate(deltaE = E3 - E1) %>% arrange(Cluster) 
  
  #### E2 - frequencies just after vaccine intro #### 
  x_imputed <- dfFNVTImputed$E1
  x_imputed <- as.matrix(round(x_imputed/sum(x_imputed), digits = 6))
  
  ##this imputation needs to see dfFVT; if both E1 and E2 are zero in NVT, then imputed from E3
  dat2_imputed <- dfImputed %>% mutate(Epoch2 = dplyr::case_when(Cluster == 16 ~ "E3",
                                                                 Cluster == 31 ~ "E3",
                                                                 TRUE  ~ "E2")) %>% 
    dplyr::select(Cluster,vaccinetype,Epoch2) %>% bind_rows(SCE2) %>% 
    arrange(Cluster) %>% left_join(SC_freq_df) %>%
    dplyr::select(btuD_2:gatA_3) 
  
  ### Get the matrix and the SC for the pre-vaccine epoch "E1"
  df_preV <- SCE1 %>% left_join(SC_freq_df)
  SC_freq_preV <- as.matrix(df_preV %>% mutate(SC_freq=SC_n/sum(SC_n)) %>% dplyr::select(SC_freq))
  SC_COG_preV  <- as.matrix(t(df_preV %>% dplyr::select(btuD_2:gatA_3)))
  
  #### Get e_l ####
  el <- SC_COG_preV %*% SC_freq_preV
  
  #### fitness function pieces #### 
  fl_imp <- t(dat2_imputed) %*% x_imputed
  
  ##for predictions
  df_postV <- data.frame(dat2_imputed)
  SC_COG_postV <- as.matrix(t(df_postV %>% dplyr::select(btuD_2:gatA_3)))
  
  SCE3 <- dfFNVTImputed %>% dplyr::select(Cluster, vaccinetype)
  SC_freq_postV_obs <- SCE3 %>% mutate(Epoch2 = "E3") %>% 
    left_join(SC_freq_df) %>%
    mutate(SC_freq=SC_n/sum(SC_n, na.rm = T)) %>% 
    dplyr::select(Cluster,vaccinetype, SC_freq) %>% 
    mutate(SC_freq = replace_na(SC_freq, 0))
  
  ## build vaccine table 
  vaccineT <- df %>%
    distinct(Cluster, vaccinetype) %>%
    subset(vaccinetype == "VT") %>%
    dplyr::rename(W = vaccinetype) %>%
    full_join(
      df %>% distinct(Cluster, vaccinetype) %>% subset(vaccinetype == "NVT"),
      by = "Cluster"
    ) %>%
    tidyr::unite("vaccine", W:vaccinetype, na.rm = TRUE)
  
  ## Predict postV frequencies
  SC_freq_postV_pred <- QP(SC_COG_postV, el)
  
  SC_freq_postV_obs <- SC_freq_postV_obs %>% 
    mutate(SC_pred = SC_freq_postV_pred) %>%
    left_join(vaccineT, by = "Cluster")
  
  ###change below to make it shown on the figures
  stats <- summary(lm(SC_freq_postV_obs$SC_freq ~SC_freq_postV_obs$SC_pred))
  ars   <- round(stats$adj.r.squared, digits = 3)
  sseE  <- round(sse(SC_freq_postV_obs$SC_freq, SC_freq_postV_obs$SC_pred), digits = 3)
  rmseE <- round(rmse(SC_freq_postV_obs$SC_freq, SC_freq_postV_obs$SC_pred), digits = 3)
  
  outlier_unw <- SC_freq_postV_obs %>% mutate(diff = abs(SC_freq - SC_pred))
  outlier_unw <- outlier_unw %>% filter(diff %in% boxplot(outlier_unw$diff, plot = FALSE)$out)
  
  plot5a <- ggplot(SC_freq_postV_obs, 
                   aes(x = SC_pred, y = SC_freq, colour = vaccine)) + 
    geom_segment(aes(x=0,xend=0.16,y=0,yend=0.16), 
                 color="black",alpha=.7,lwd=0.5,lty=3) +
    theme(legend.position = "none") + theme_classic() + 
    geom_smooth(method='lm', color="#899DA4" ,
                formula=y~x, alpha=0.3, lwd=.6, 
                fullrange=T, linetype="blank", show.legend=F) +
    annotate(geom = "text", x=0.145, y =0.15, 
             label = "1:1 line", angle = 45, size = 3) + 
    geom_point(size=3, alpha = 0.9) +  
    scale_x_continuous("Predicted Prevalence (NFDS)")+
    scale_y_continuous("Observed Prevalence") +
    coord_fixed(ratio = 1, xlim=c(0,0.16), ylim=c(0,0.16)) + 
    scale_colour_manual(values = c("#143c77","darkorchid4","darkred")) + 
    theme(legend.position = "none", 
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(size = 15))  + 
    annotate("text", x=0.001, y=0.135, size=5,hjust = 0,
             label=paste("SSE = ", sseE, "\nRMSE = ", 
                         rmseE, "\nAdj. R2 = ", ars)) + 
    geom_text_repel(aes(label = paste("SC", Cluster, sep = "-")), data = outlier_unw, size = 4)


  #############
  # PREDICTION CODE: WEIGHTED 
  #############
  #apply weights 
  dfFVT_w <- df4 %>% dplyr::select(Cluster,vaccinetype, Epoch2, weightadjust) %>% 
    group_by(Epoch2,Cluster, vaccinetype) %>% summarise(cnt = sum(weightadjust)) %>%ungroup() %>%
    group_by(Epoch2)%>%
    mutate(freq = round(cnt/sum(cnt), digits = 6)) %>% ungroup()%>%
    dplyr::select(Epoch2, Cluster,vaccinetype, freq) %>% 
    spread(Epoch2, freq, fill = 0) %>% arrange(Cluster,vaccinetype)
  
  SC_freq_denom <- df4 %>% dplyr::select(Cluster, vaccinetype, Epoch2,weightadjust, 
                                         btuD_2:gatA_3) %>%
    arrange(Cluster) %>% group_by(Cluster,vaccinetype,Epoch2) %>%
    summarize(denominator=sum(weightadjust))
  
  group_names <- colnames(df4)
  group_names <- group_names[!(group_names %in% c("Cluster","vaccinetype","Epoch2"))]
  
  SC_freq_df_w <- df4 %>% dplyr::select(Cluster, vaccinetype, Epoch2,weightadjust, 
                                        btuD_2:gatA_3) %>%
    arrange(Cluster) %>% dplyr::select(-weightadjust) %>% pivot_longer(-c(Cluster,vaccinetype,Epoch2)) %>%
    group_by(Cluster,vaccinetype,Epoch2, name) %>%
    summarize(sum_value=sum(value)) %>%
    left_join(SC_freq_denom) %>%
    mutate(y=sum_value/denominator)
  SC_freq_df_w$name <- factor(SC_freq_df_w$name, levels=group_names)
  
  SC_freq_df_w <- SC_freq_df_w %>% arrange(name) %>% pivot_wider(id_cols=c(Cluster,vaccinetype,Epoch2),names_from=name,values_from=y)
  SC_freq_df_w <- SC_freq_denom %>% left_join(SC_freq_df_w)
  
  #### Present at E1 ####
  SCE1_w <- dfFVT_w %>% subset(E1 > 0) %>% 
    dplyr::select(Cluster,vaccinetype) %>% 
    mutate(Epoch2 = "E1")
  
  #### NVTs Present at E1 ####
  SCE2_w <- SCE1_w %>% subset(vaccinetype == "NVT")
  
  dfImputed_w <- dfFVT_w %>% subset(vaccinetype == "NVT" & (E1 == 0 )) %>%
    dplyr::select(Cluster,vaccinetype) %>% mutate(Epoch2 = "E1", cnt=1) %>%
    dplyr::select(Epoch2,Cluster,vaccinetype,cnt) 
  
  dfFVTImputed_w <- df4 %>% dplyr::select(Cluster,vaccinetype, Epoch2, weightadjust) %>% 
    group_by(Epoch2,Cluster, vaccinetype) %>% summarise(cnt = sum(weightadjust)) %>%ungroup() %>% 
    bind_rows(dfImputed_w) %>% 
    group_by(Epoch2)%>%
    mutate(freq = round(cnt/sum(cnt), digits = 6)) %>% ungroup()%>%
    dplyr::select(Epoch2, Cluster,vaccinetype, freq) %>% 
    spread(Epoch2, freq, fill = 0) %>% arrange(Cluster,vaccinetype)
  
  dfFNVTImputed_w <- dfFVTImputed_w %>% 
    subset(vaccinetype == "NVT") %>% 
    mutate(deltaE = E3 - E1) %>% arrange(Cluster)
  
  x_imputed_w <- dfFNVTImputed_w$E1
  x_imputed_w <- as.matrix(round(x_imputed_w/sum(x_imputed_w), digits = 6))
  
  dat2_imputed_w <- dfImputed_w %>% mutate(Epoch2 = dplyr::case_when(Cluster == 31 ~ "E3", Cluster == 16 ~ "E3", TRUE ~ "E2")) %>%
    dplyr::select(Cluster,vaccinetype,Epoch2) %>% bind_rows(SCE2_w) %>% 
    arrange(Cluster) %>% left_join(SC_freq_df_w) %>%
    dplyr::select(btuD_2:gatA_3) 
  
  df_preV_w <- SCE1_w %>% left_join(SC_freq_df_w)
  SC_freq_preV_w <- as.matrix(df_preV_w %>% mutate(SC_freq= denominator/sum(denominator)) %>% dplyr::select(SC_freq))
  SC_COG_preV_w <- as.matrix(t(df_preV_w %>% dplyr::select(btuD_2:gatA_3)))
  el_w <- SC_COG_preV_w %*% SC_freq_preV_w
  
  df_postV_w <- dat2_imputed_w
  SC_COG_postV_w <- as.matrix(t(df_postV_w %>% dplyr::select(btuD_2:gatA_3)))
  
  fl_imp_w <- t(dat2_imputed_w) %*% x_imputed_w
  
  SCE3_w <- dfFNVTImputed_w %>% dplyr::select(Cluster, vaccinetype)
  SC_freq_postV_obs_w <- SCE3_w %>% mutate(Epoch2 = "E3") %>% 
    left_join(SC_freq_df_w) %>%
    mutate(SC_freq=denominator/sum(denominator, na.rm = T)) %>% 
    dplyr::select(Cluster,vaccinetype, SC_freq) %>% 
    mutate(SC_freq = replace_na(SC_freq, 0))
  
  ## Predict postV frequencies
  SC_freq_postV_pred_w <- QP(SC_COG_postV_w, el_w)
  
  SC_freq_postV_obs_w <- SC_freq_postV_obs_w %>% 
    mutate(SC_pred = SC_freq_postV_pred_w) %>%
    left_join(vaccineT, by = "Cluster")
  
  ###change below to make it shown on the figures
  stats_w <- summary(lm(SC_freq_postV_obs_w$SC_freq ~SC_freq_postV_obs_w$SC_pred))
  ars_w   <- round(stats_w$adj.r.squared, digits = 3)
  sseE_w  <- round(sse(SC_freq_postV_obs_w$SC_freq, SC_freq_postV_obs_w$SC_pred), digits = 3)
  rmseE_w <- round(rmse(SC_freq_postV_obs_w$SC_freq, SC_freq_postV_obs_w$SC_pred), digits = 3)
  
  outlier_w <- SC_freq_postV_obs_w %>% mutate(diff = abs(SC_freq - SC_pred))
  outlier_w <- outlier_w %>% filter(diff %in% boxplot(outlier_w$diff, plot = FALSE)$out)
  
  plot5b <- ggplot(SC_freq_postV_obs_w, 
                   aes(x = SC_pred, y = SC_freq, colour = vaccine)) + 
    geom_segment(aes(x=0,xend=0.16,y=0,yend=0.16), 
                 color="black",alpha=.7,lwd=0.5,lty=3) +
    theme(legend.position = "none") + theme_classic() + 
    geom_smooth(method='lm', color="#899DA4" ,
                formula=y~x, alpha=0.3, lwd=.6, 
                fullrange=T, linetype="blank", show.legend=F) +
    annotate(geom = "text", x=0.145, y =0.15, 
             label = "1:1 line", angle = 45, size = 3) + 
    geom_point(size=3, alpha = 0.9) +  
    scale_x_continuous("Predicted Prevalence (NFDS)")+
    scale_y_continuous("Observed Prevalence") +
    coord_fixed(ratio = 1, xlim=c(0,0.16), ylim=c(0,0.16)) + 
    scale_colour_manual(values = c("#143c77","darkorchid4","darkred")) + 
    theme(legend.position = "none", 
          axis.text = element_text(colour = "black", size = 10),axis.title = element_text(size = 15)) + 
    annotate("text", x=0.001, y=0.135, size=5,hjust = 0,
             label=paste("SSE = ", sseE_w, "\nRMSE = ", 
                         rmseE_w, "\nAdj. R2 = ", ars_w)) + 
    geom_text_repel(aes(label = paste("SC", Cluster, sep = "-")), data = outlier_w, size = 4)
  
  # combine and save
  fig5 <- ggarrange(plot5a, plot5b, ncol = 2)
  ggsave(paste0("Fig5_", suffix, ".pdf"), fig5, width = 12, height = 6)
  
  return(
    list(
      metrics = tibble(
        scenario   = suffix,
        adj_r2_unw = ars,
        sse_unw    = sseE,
        rmse_unw   = rmseE,
        adj_r2_w   = ars_w,
        sse_w      = sseE_w,
        rmse_w     = rmseE_w
      ),
      plot_unw = plot5a,
      plot_w   = plot5b
    )
  )
}

###############################################################################
# 5. RUN SCENARIOS
###############################################################################

# main
res_main   <- run_prediction_sensitivity(df2_master, df4_master,
                                         gpsc_cutoff = 5,
                                         sero3_as_nvt = FALSE,
                                         suffix = "gpsc5")

# GPSC > 3 --> Supplemental Fig S4
res_gpsc3  <- run_prediction_sensitivity(df2_master, df4_master,
                                         gpsc_cutoff = 3,
                                         sero3_as_nvt = FALSE,
                                         suffix = "supplemental_gpsc3")

# GPSC > 10 --> Supplemental Fig S5
res_gpsc10 <- run_prediction_sensitivity(df2_master, df4_master,
                                         gpsc_cutoff = 10,
                                         sero3_as_nvt = FALSE,
                                         suffix = "supplemental_gpsc10")

# serotype 3 → NVT in vaccinetype --> Supplemental Fig S6
res_s3nvt  <- run_prediction_sensitivity(df2_master, df4_master,
                                         gpsc_cutoff = 5,
                                         sero3_as_nvt = TRUE,
                                         suffix = "supplemental_gpsc5_sero3NVT")
