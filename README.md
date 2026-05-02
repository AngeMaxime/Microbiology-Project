# 🦠 Task Code 2.1 — Microbiology: Bacterial Growth Curve Analysis

![R](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r&logoColor=white)
![RMarkdown](https://img.shields.io/badge/Report-RMarkdown-blue?style=flat)
![HackBio](https://img.shields.io/badge/HackBio-Internship%202025-orange?style=flat)
![License](https://img.shields.io/badge/License-MIT-green?style=flat)

---

## 🧭 Introduction

This project is part of the **HackBio Internship 2025 — Task Code 2.1: Microbiology**.

In microbiology, bacterial growth follows a well-known sigmoidal pattern — passing through a **lag phase**, an **exponential growth phase**, and finally a **stationary phase** where the population stabilises at its **carrying capacity (K)**. The time it takes a bacterial strain to reach this carrying capacity is a key indicator of its growth performance and competitive fitness.

This analysis uses the **MCGC (Microbial Growth Curve Collection)** dataset, which contains **OD600 optical density readings** measured over time for bacterial strains grown in a 96-well plate. Each strain has been tested in two genetic configurations:

- 🔴 **Knock-OUT (-)** — the gene of interest has been deleted or disabled
- 🔵 **Knock-IN (+)** — the gene has been inserted or activated

| 📁 Resource | 🔗 Link |
|---|---|
| 📊 Dataset | [mcgc.tsv](https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/mcgc.tsv) |
| 📖 Dataset Description | [mcgc_METADATA.txt](https://github.com/HackBio-Internship/2025_project_collection/blob/main/Python/Dataset/mcgc_METADATA.txt) |

---

## 🎯 Purpose

The purpose of this analysis is to compare the **growth dynamics** of Knock-OUT and Knock-IN bacterial strains and determine whether the presence or absence of the gene of interest has a statistically significant effect on the time required to reach carrying capacity.

---

## 🔍 What We Want to Do

- 📥 Load and inspect the MCGC dataset and map each well to its strain and mutant type
- 🔄 Reshape the data into long format suitable for visualisation
- 📈 Plot growth curves (OD600 vs Time) for each strain with KO and KI overlaid
- ⚙️ Define a function to compute the **time to reach carrying capacity** for each well
- 🔵 Generate a **scatter plot** comparing time to carrying capacity between KO and KI
- 📦 Generate a **box plot** comparing the distribution of time to carrying capacity between KO and KI
- 📐 Run a **Welch two-sample t-test** to determine whether the difference is statistically significant
- 💬 Interpret all observations as inline comments in the code

---

## 🔬 Methods

### 📂 1. Data Import & Plate Layout
The dataset is loaded directly from GitHub as a tab-separated file. Each row is a time point (measured every 15 minutes), and each column (A1–A12, B1–B12, C1–C12) represents one well in a 96-well plate.

A **metadata table** is manually built to map each well to its strain and mutant type based on the plate layout:

| 🧫 Plate Row | Wells 1–6 | Wells 7–12 |
|---|---|---|
| **A** | Strain A — 🔴 Knock-OUT (-) | Strain A — 🔵 Knock-IN (+) |
| **B** | Strain B — 🔴 Knock-OUT (-) | Strain B — 🔵 Knock-IN (+) |
| **C** | Strain C — 🔴 Knock-OUT (-) | Strain C — 🔵 Knock-IN (+) |

💾 Three files are exported at this stage: `Raw_data_microbiology.csv`, `Meta_data.csv`, and `Long_format.csv`.

### 🔄 2. Data Reshaping
The wide-format dataset (one column per well) is pivoted into **long format** using `pivot_longer()`, then joined with the metadata table. Each row in the resulting dataframe represents one OD600 reading at a specific time point for a specific well.

### 📈 3. Growth Curve Visualisation
For each strain, the OD600 vs Time growth curves of KO and KI replicates are **overlaid** on the same panel:
- 〰️ Thin transparent lines = individual well replicates
- ➖ Bold solid lines = group mean per mutant type (computed with `stat_summary()`)
- 🗂️ Faceted by strain (one panel per strain)

💾 Saved as `growth_plot.png`.

### ⚙️ 4. Time to Carrying Capacity Function
A custom function `time_to_carrying_capacity()` is defined:

```r
time_to_carrying_capacity <- function(time_vec, od_vec, threshold = 0.95) {
  K   <- max(od_vec, na.rm = TRUE)          # carrying capacity = max OD600
  idx <- which(od_vec >= threshold * K)[1]  # first time >= 95% of K
  if (is.na(idx)) return(NA)
  return(time_vec[idx])
}
```

Using a **95% threshold** (instead of the exact maximum) avoids noise sensitivity at the plateau. The function is applied to every well using `group_by()` + `summarise()`.

### 🔵 5. Scatter Plot
A scatter plot is generated with:
- X axis = Mutant type (KO vs KI)
- Y axis = Time to carrying capacity (minutes)
- Points colored by mutant type, shaped by strain
- A **crossbar** showing the group mean

💾 Saved as `Scatter.png`.

### 📦 6. Box Plot
A box plot is generated showing the distribution of time to carrying capacity:
- Fill color = mutant type (KO vs KI)
- Jittered points colored by strain (A, B, C)
- Box = interquartile range (IQR), line = median

💾 Saved as `boxplot.png`.

### 📐 7. Statistical Test
A **Welch two-sample t-test** is applied (does not assume equal variances):
- **H₀:** Mean time to carrying capacity is the **same** for KO and KI
- **H₁:** Mean time to carrying capacity **differs** between KO and KI

The p-value is automatically interpreted: if p < 0.05, the difference is statistically significant. An annotated boxplot with significance stars (`*`, `**`, `***`, `ns`) is produced using `ggpubr::stat_compare_means()` and saved as `boxplot_pval.png`.

---

## 📊 Results

### 📈 Growth Curves
All three strains (A, B, C) exhibit a **classic sigmoidal growth pattern**:

| 🕐 Phase | Description |
|---|---|
| **1️⃣ Lag phase** | Flat OD600 at the start — bacteria adapting to the environment |
| **2️⃣ Log phase** | Rapid exponential increase in OD600 — active cell division |
| **3️⃣ Stationary phase** | OD600 plateaus at the carrying capacity (K) — growth stabilises |

🔵 Knock-IN (+) strains generally reach carrying capacity **faster** and/or at a **higher OD600** than 🔴 Knock-OUT (-) strains, suggesting the inserted gene confers a growth advantage.

### ⏱️ Time to Carrying Capacity
- 🔵 **Knock-IN (+)** points cluster at **lower** time values → KI strains reach their plateau faster
- 🔴 **Knock-OUT (-)** points cluster at **higher** time values → KO strains take longer
- 📏 The **IQR of KO is wider**, indicating greater variability in growth timing when the gene is absent — the gene also contributes to growth **consistency**

### 📐 Statistical Test
The Welch t-test determines whether the observed difference is statistically significant (p < 0.05). A significant result indicates that the gene of interest plays a **biologically meaningful role** in accelerating bacterial growth — likely encoding a function related to **metabolism, nutrient uptake, or cell division**.

---

## 💡 What We Can Retain

> 🔵 **Knock-IN (+) strains consistently reach carrying capacity faster** than Knock-OUT (-) strains across all three strains, confirming that the inserted gene confers a meaningful growth advantage.

> 📏 **The gene also contributes to growth consistency** — its absence (KO) leads to more variable and unpredictable growth timing, as shown by the wider IQR in the box plots.

> 📐 **The statistical test (t-test)** confirms whether this difference is significant. A p-value < 0.05 means we can confidently conclude that the gene's presence or absence has a real biological impact on bacterial population dynamics, not just random variation.

> 🧬 **Biological implication:** The knocked-in gene likely encodes a protein involved in a growth-promoting pathway — its presence accelerates the rate at which bacteria multiply and fill their environment to maximum capacity.

---

## 🧠 Skills Learned

- 📂 Loading and inspecting `.tsv` files in R with `read.delim()`
- 🗂️ Manually building metadata tables to annotate experimental layouts
- 🔄 Reshaping data from wide to long format with `tidyr::pivot_longer()`
- 🔗 Joining datasets with `dplyr::left_join()`
- 📈 Plotting overlaid growth curves with `ggplot2` using `stat_summary()` for mean lines
- ⚙️ Writing custom R functions with a biological threshold parameter
- 🧮 Applying grouped summaries with `dplyr::group_by()` + `summarise()`
- 🎨 Visualising group comparisons with scatter plots and box plots in `ggplot2`
- 📐 Running and interpreting Welch two-sample t-tests in R
- ⭐ Annotating plots with statistical significance using `ggpubr::stat_compare_means()`
- 💾 Exporting plots as high-resolution PNG files with `ggsave()`

---

## 👤 Authors

**Ange Maxime TCHOUTANG**
🎓 HackBio Internship 2025 — Task Code 2.1: Microbiology

---

## 📦 Dependencies

| 📦 Package | 🎯 Purpose |
|---|---|
| `tidyverse` | Data manipulation (`dplyr`, `tidyr`) + `ggplot2` |
| `ggplot2` | All growth curve and comparative visualisations |
| `patchwork` | Combining multiple plots into panels |
| `scales` | Axis label formatting helpers |
| `ggpubr` | Statistical significance annotation on plots |

---

## 🔗 Reference

*HackBio Internship 2025 — https://github.com/HackBio-Internship/2025_project_collection*
