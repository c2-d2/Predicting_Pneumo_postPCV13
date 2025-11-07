# 🧬 Predicting Pneumo post-PCV13

This repository contains the data and analysis scripts used in the study:  
**"Predicting pneumococcal population composition after vaccine introduction"**  
Published in *mBio* (2024)  
🔗 [https://journals.asm.org/doi/full/10.1128/mbio.03355-23](https://journals.asm.org/doi/full/10.1128/mbio.03355-23)

---

## 📘 Overview

This project investigates how pneumococcal population composition changes following pneumococcal conjugate vaccine (PCV13) introduction, using genomic data and a **negative frequency-dependent selection (NFDS)** model.  
The scripts implement data preparation, weighting, and model-based prediction of strain (GPSC) frequencies.

---

## 📂 Repository Structure

```
Predicting_Pneumo_postPCV13/
├── data/                     # All input data files (metadata, weights, GPSC cluster assignments, etc.)
├── output/                   # All output files and figures
├── 01_accessory_filter.R     # Step 1: Filter and convert Roary accessory genome output
├── 02_weight_eval.R          # Step 2: Generate serotype-specific weights and evaluate correlation
├── 03_prediction_model.R     # Step 3: Apply NFDS model to predict post-PCV13 strain frequencies
└── README.md                 # Project documentation
```

---

## 🧭 How to Use

### 1️⃣ Clone or Download this Repository

```bash
git clone https://github.com/c2-d2/Predicting_Pneumo_postPCV13.git
```

Or click **Code → Download ZIP** on GitHub and extract the folder.

---

### 2️⃣ Open the Project in RStudio

1. Create a new **R Project** using this folder as the working directory.  
2. Ensure your working directory points to the root of the downloaded repo.

---

### 3️⃣ Install Required R Packages

You can install all dependencies at once:

```r
install.packages(c("tidyverse", "readxl", "ggpubr", "cowplot", "quadprog", "Metrics", "EnvStats"))
```

---

### 4️⃣ Run the Scripts in Order

Run the following scripts sequentially:

1. `presence-absence-matrix.R`  
2. `Correlation_analysis.R`  
3. `NFDS_prediction_function_sensitivity.R`

Each script generates intermediate outputs required for the next step.  
All data are automatically loaded from the `data/` folder.

---

### 5️⃣ Outputs

1. Filtered accessory gene tables  
2. Serotype-specific weight figures  
3. Model predictions and evaluation plots (e.g., `Fig4.pdf`, `Fig5.pdf`)

---

## 🧩 Data

All datasets required to reproduce the analysis are included in the `/data` folder:

- `CDC98_18_IPD_UScarriage-accessory-gene-with5-95percent.csv`  
- `cdc-weight-final-022022.xlsx`  
- `cdcIPD-GPSC-assigned_clusters2.csv`  
- `cdc-all-metadata.xlsx`  
- `cdc98-18_IPD_UScarriage_allmetadata.xlsx`  

**Note:**  
For Script 1 input (`gene_presence_absence.csv`), the file is too large for GitHub and is hosted externally on Dropbox.  
You can download it here:  
🔗 [Dropbox Link](https://www.dropbox.com/scl/fo/8qj9mjslhagh6xsizupeh/AEw8UVQgkHURlqoVLTfEkbs?rlkey=prshrezrarnz3iaoqyqxt9bwl&st=vlmmdmg3&dl=0)

---

## 🧠 Notes

- Scripts were developed and tested in **R 4.2.2**.  
- Figures and outputs can be reproduced exactly if the directory structure is preserved.  
- Sensitivity analyses can be conducted by modifying:
  - **GPSC cluster size cutoffs** (e.g., >3, >5, >10 isolates)  
  - **Serotype 3 classification**, optionally treating it as a non-vaccine type (NVT)

---

## 👩‍🔬 Author

**Xueting Qiu, PhD**  
Harvard T.H. Chan School of Public Health (former postdoc)  

📧 xuetingqiu@hsph.harvard.edu  
📧 qiuxueting1402@gmail.com  

---

## 📜 Citation

If you use this code or data, please cite:

> Qiu X., et al. (2023). *Predicting pneumococcal population composition after vaccine introduction.*  
> *mBio* 14(2): e03355-23.  
> DOI: [10.1128/mbio.03355-23](https://doi.org/10.1128/mbio.03355-23)
