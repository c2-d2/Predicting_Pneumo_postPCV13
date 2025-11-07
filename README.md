# Predicting Pneumo post-PCV13

This repository contains the data and analysis scripts used in the study:  
**"Predicting pneumococcal population composition after vaccine introduction"**  
Published in *mBio* (2023)  
🔗 [https://journals.asm.org/doi/full/10.1128/mbio.03355-23](https://journals.asm.org/doi/full/10.1128/mbio.03355-23)

---

## 📘 Overview
This project investigates how pneumococcal population composition changes following pneumococcal conjugate vaccine (PCV13) introduction, using genomic data and a negative frequency-dependent selection (NFDS) model.  
The scripts implement data preparation, weighting, and model-based prediction of strain (GPSC) frequencies.

---

## 📂 Repository Structure

Predicting_Pneumo_postPCV13/
├── data/ # All input data files (metadata, weights, GPSC cluster assignments, etc.)
├── output/ # All output files and figures
├── 01_accessory_filter.R # Step 1: Filter and convert Roary accessory genome output
├── 02_weight_eval.R # Step 2: Generate serotype-specific weights and evaluate correlation
├── 03_prediction_model.R # Step 3: Apply NFDS model to predict post-PCV13 strain frequencies
└── README.md # Project documentation

---

## 🧭 How to Use

1. **Clone or Download this repository**
   ```bash
   git clone https://github.com/c2-d2/Predicting_Pneumo_postPCV13.git
or click Code → Download ZIP on GitHub and extract the folder.

2. Open the project in RStudio

Create a new R Project using this folder as the working directory.

Ensure your working directory points to the root of the downloaded repo.

Install required R packages
(You can install them all at once)

install.packages(c("tidyverse", "readxl", "ggpubr", "cowplot", "quadprog", "Metrics", "EnvStats"))

3. Run the scripts in order

    1.presence-absence-matrix.R
    2.Correlation_analysis.R
    3.NFDS prediction_function_sensitivity.R

Each script is designed to generate intermediate results for the next one.
All data are automatically loaded from the data/ folder.

4. Outputs

    1) Filtered accessory gene tables
    
    2) Serotype-specific weight figures
    
    3) Model predictions and evaluation plots (e.g., Fig4.pdf, Fig5.pdf)

🧩 Data
All datasets required to reproduce the analysis are included in the /data folder:

      CDC98_18_IPD_UScarriage-accessory-gene-with5-95percent.csv
      
      cdc-weight-final-022022.xlsx
      
      cdcIPD-GPSC-assigned_clusters2.csv
      
      cdc-all-metadata.xlsx
      
      cdc98-18_IPD_UScarriage_allmetadata.xlsx
      
      For the input of R scripts 1 -- gene_presence_absence.csv, the file is too big, so it is hosted in dropbox. You can download it here: https://www.dropbox.com/scl/fo/8qj9mjslhagh6xsizupeh/AEw8UVQgkHURlqoVLTfEkbs?rlkey=prshrezrarnz3iaoqyqxt9bwl&st=vlmmdmg3&dl=0


🧠 Notes
The scripts were developed and tested in R 4.2.2.

Figures and outputs can be reproduced exactly if directory structure is preserved.

Sensitivity analyses can be conducted by modifying GPSC cluster size cutoffs and vaccine type definitions in the final prediction script.

👩‍🔬 Author
Xueting Qiu, PhD, MD-equivalent
Senior Scientist, Ginkgo Biosecurity
Harvard T.H. Chan School of Public Health (former postdoc)
📧 xuetingqiu@hsph.harvard.edu

📜 Citation
If you use this code or data, please cite:

Qiu X., et al. (2023). Predicting pneumococcal population composition after vaccine introduction. mBio 14(2):e03355-23.
DOI: 10.1128/mbio.03355-23
