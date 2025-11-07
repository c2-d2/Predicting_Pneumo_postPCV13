###############################################################################
# Accessory Genome Filtering and 0/1 Matrix Converting Analysis
#
# Description:
#   This script processes the gene presence/absence output from Roary pangenome analysis to:
#     1. Filter COGs based on frequency and length
#     2. Convert presence/absence data into a binary (0/1) matrix
#     3. Join with metadata tables (IPD and US carriage)
#     4. Export accessory-genome table with metadata and 0/1 matrix
#
# Notes:
#   - Adjust 'roary_start_col' depending on Roary version:
#       * Old Roary output:  start column = 12
#       * New Roary output:  start column = 15
#
# Author: Xueting Qiu
# Original:   2019-2020 
# Refactored: Nov.2025
###############################################################################

library(tidyverse)
library(readxl)
library(writexl)
library(scatterD3)
library(here)

# ---- params ---------------------------------------------------------------
roary_start_col <- 15      # 15 for new Roary; 12 for old Roary
trim_percent   <- 0.10     # trimming low-frequency COGs
trim_size_nt   <- 150      # min COG length
upper_prop     <- 0.95     # upper bound for accessory definition
lower_prop     <- 0.05     # lower bound for accessory definition
n_meta_cols    <- 8        # how many columns in metadata to keep in front

# ---- read data ------------------------------------------------------------
gene_presence_absence <- readr::read_csv(here("data/gene_presence_absence.csv"))
cdcUScarriage_meta    <- readxl::read_excel(here("data/cdc98-18_IPD_UScarriage_allmetadata.xlsx"))
#metadata has the accesssion numbers for all the sequences. 

# ---- basic stats ----------------------------------------------------------
unfiltered_pangenome_size <- nrow(gene_presence_absence) - 1
taxa.number <- ncol(gene_presence_absence) - (roary_start_col - 1)

print(unfiltered_pangenome_size)
print(taxa.number)

# ---- diagnostics: unfiltered ----------------------------------------------
hist(gene_presence_absence$`No. isolates`,
     main = "COG Frequencies (unfiltered)",
     xlab  = "Number of Taxa")
hist(gene_presence_absence$`Avg group size nuc`,
     main = "Distribution of COG sizes (nt, unfiltered)",
     xlab  = "COG size (nt)")

ggplot(gene_presence_absence, aes(`No. isolates`)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "No. of Isolates containing COG", y = "Frequency")

# ---- filter by length + frequency (Roary-level) ---------------------------
df <- gene_presence_absence %>%
  filter(!(`Avg group size nuc` < trim_size_nt &
             `No. isolates` < max(`No. isolates`) * trim_percent))

genename <- df[1]  # first column is the gene/group name

# diagnostics on filtered
hist(df$`No. isolates`,
     main = "COG Frequencies (filtered)",
     xlab  = "Number of Taxa")
hist(df$`Avg group size nuc`,
     main = "Distribution of COG sizes (nt, filtered)",
     xlab  = "COG size (nt)")

scatterD3(df$`Avg group size nuc`,
          df$`No. isolates`,
          xlab = "Average Length (nt)",
          ylab = "Number of Taxa",
          axes_font_size = "90%")

# zoom small
df2 <- df[df$`Avg group size nuc` < 500, ]
scatterD3(df2$`Avg group size nuc`,
          df2$`No. isolates`,
          xlab = "Average Length (nt)",
          ylab = "Number of Taxa",
          axes_font_size = "90%")
df2 <- NULL

# ---- build presence/absence 0/1 matrix (shared for all metadata) ----------
# keep only isolate columns
pa <- df[, roary_start_col:ncol(df)]

# convert to 0/1
pa[] <- lapply(pa, as.character)
pa[] <- lapply(pa, function(x) ifelse(is.na(x) | x == "", 0L, 1L))

# add gene names, then transpose so isolates are rows
df_col <- cbind(genename, pa)
df_t <- as.data.frame(t(df_col), stringsAsFactors = FALSE)
colnames(df_t) <- df_t[1, ]
df_t <- df_t[-1, ]
df_t <- tibble::rownames_to_column(df_t, var = "Accession")
df_t[] <- lapply(df_t, type.convert, as.is = TRUE)

# ---- helper function to create accessory table for a given metadata -------
make_accessory_table <- function(meta_df,
                                 df_t,
                                 id_col_in_meta = "seq_filename",  # change if your metadata uses another name
                                 keep_meta_cols = n_meta_cols,
                                 lower = lower_prop,
                                 upper = upper_prop,
                                 out_path = NULL,
                                 plot_hist = TRUE) {
  
  # 1) standardize accession name in metadata
  meta_df <- meta_df %>% rename(Accession = !!id_col_in_meta)
  
  # 2) join metadata with presence/absence matrix
  comb <- meta_df %>%
    inner_join(df_t, by = "Accession")
  
  # 3) pull out PA-only part (after metadata columns)
  pa_with_meta <- comb[, (keep_meta_cols + 1):ncol(comb)]
  
  # 4) calculate frequencies and filter accessory genes
  taxa.number <- nrow(pa_with_meta)
  col_freq    <- colSums(pa_with_meta)
  
  df_acc <- pa_with_meta[
    ,
    col_freq < round(upper * taxa.number, 0) &
      col_freq > round(lower * taxa.number, 0)
  ]
  
  # 5) optional histogram check on COG frequency distribution
  if (isTRUE(plot_hist)) {
    hist(colSums(df_acc) / taxa.number,
         main = "Accessory Genome COG Frequencies",
         xlab  = "COG frequency")
  }
  
  # 6) join metadata back on the left
  out <- cbind(comb[, 1:keep_meta_cols], df_acc)
  
  # 7) write file if asked
  if (!is.null(out_path)) {
    readr::write_csv(out, out_path)
  }
  
  return(out)
}


# ---- run for CDC US carriage metadata -------------------------------------
out_ipd_carriage <- make_accessory_table(
  meta_df        = cdcUScarriage_meta,
  df_t           = df_t,
  id_col_in_meta = "seq_filename",   # or "accession" depending on your sheet
  keep_meta_cols = 8,
  lower          = lower_prop,
  upper          = upper_prop,
  out_path       = here("data/CDC98_18_IPD_UScarriage-accessory-gene-with5-95percent.csv")
)

##### Notes:
# CDC IPD + US carriage sequence n = 14046; COGs n=3218; COGs columns: (butD_2:gatA_3)
#check last column name: dplyr::last(colnames(out_ipd_carriage))