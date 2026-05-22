# MATLAB-Based Numerical Simulation of BCR-ABL1 Expression Dynamics and Targeted Therapy Response in Chronic Myeloid Leukemia

[![MATLAB](https://img.shields.io/badge/MATLAB-Numerical%20Modeling-orange?style=flat-square)](../MATLAB_Code/CML_Numerical_Project.m)
[![Cancer Modeling](https://img.shields.io/badge/Cancer%20Modeling-CML-red?style=flat-square)](#)
[![Computational Biology](https://img.shields.io/badge/Computational%20Biology-Gene%20Expression-blue?style=flat-square)](#)

## Project Links

- **Main MATLAB Code:** [CML_Numerical_Project.m](../MATLAB_Code/CML_Numerical_Project.m)
- **Portfolio README:** [Back to main profile](../README.md)

---

## Project Overview

This project presents a MATLAB-based computational and numerical simulation of **BCR-ABL1 expression dynamics** in **Chronic Myeloid Leukemia (CML)**. It combines biological interpretation, gene expression analysis, and numerical methods to model how leukemic activity changes under targeted therapy.

The project focuses on the biological role of **BCR-ABL1**, the oncogenic fusion driver of CML. BCR-ABL1 produces abnormal tyrosine kinase activity that supports uncontrolled leukemic cell growth and survival. The simulation represents how this disease-related activity may decrease under a tyrosine kinase inhibitor-like treatment effect.

The project also integrates transcriptomic analysis inspired by the public GEO dataset **GSE33075**, which includes healthy, pre-treatment CML, and post-imatinib treatment samples.

---

## Project Aim

The aim of this project was to build a simple but biologically meaningful MATLAB workflow that can:

- Simulate CML-related molecular and cellular dynamics over time.
- Compare normal, untreated CML, and treated CML behavior.
- Apply numerical methods to solve ordinary differential equations.
- Analyze gene expression differences between disease and reference groups.
- Visualize treatment-response behavior using scientific plots.

---

## Biological Background

Chronic Myeloid Leukemia is strongly associated with the **Philadelphia chromosome**, which creates the **BCR-ABL1 fusion gene**. This fusion gene produces an abnormal tyrosine kinase that activates signaling pathways involved in cell proliferation and survival.

Targeted therapies such as imatinib are designed to inhibit BCR-ABL1 kinase activity. In a simplified computational model, treatment can therefore be represented as a reduction in leukemic cell survival or molecular activity over time.

---

## Mathematical Model

The MATLAB script includes a simplified CML model using three state variables:

- **x(t):** leukemic cell population
- **y(t):** healthy cell population
- **m(t):** BCR-ABL1-related mRNA expression signal

The model describes:

- leukemic cell growth and death,
- healthy cell growth,
- competition between leukemic and healthy cells,
- mRNA production related to leukemic burden,
- mRNA degradation,
- treatment-induced reduction in leukemic cell growth/survival.

The treatment parameter is represented as **u**, where larger values indicate stronger treatment effect.

---

## Numerical Methods Used

The project applies three numerical approaches:

| Method | Purpose |
|---|---|
| **Euler’s Method** | Step-by-step approximation of the ODE system |
| **RK4 Method** | More accurate fixed-step numerical simulation |
| **MATLAB ode45** | Built-in adaptive ODE solver used as a reference |

The comparison between Euler, RK4, and ode45 demonstrates how different numerical methods approximate the same biological system.

---

## Transcriptomic Data Analysis Component

The code includes a gene expression analysis workflow that can read a CSV or Excel expression matrix. If no dataset is found, the script automatically creates a small demo dataset so the project can run without errors.

The analysis includes:

- expression data loading,
- metadata creation,
- missing-value handling,
- optional log2 transformation,
- differential expression analysis,
- Benjamini-Hochberg adjusted p-values,
- classification of genes as upregulated, downregulated, or not significant,
- output of DEG results as an Excel file.

---

## Main Outputs

When the MATLAB script is run, it creates a folder called `Project_Output` containing tables, summaries, and figures.

### Generated Result Files

| Output | Description |
|---|---|
| `Differential_Expression_Results.xlsx` | Gene-level logFC, p-value, adjusted p-value, and regulation status |
| `Treatment_Scenario_Summary.xlsx` | Final model values under different treatment strengths |
| `Project_Summary.txt` | Short biological and numerical summary of the analysis |

### Generated Figures

| Figure | Description |
|---|---|
| `Figure_1_sample_expression_boxplot.png` | Expression distribution across samples |
| `Figure_2_sample_mean_expression.png` | Mean expression per sample |
| `Figure_3_volcano_plot.png` | Differential expression volcano plot |
| `Figure_4_top_DE_genes_barplot.png` | Top differentially expressed genes |
| `Figure_5_top_genes_heatmap.png` | Heatmap of top genes |
| `Figure_6_CML_dynamics_RK4.png` | CML model dynamics using RK4 |
| `Figure_7_numerical_methods_comparison.png` | Euler, RK4, and ode45 comparison |
| `Figure_8_treatment_scenarios.png` | Effect of treatment strength on leukemic cells |

---

## How to Run the Project

1. Download or open the MATLAB script:

```text
MATLAB_Code/CML_Numerical_Project.m
```

2. Open the file in MATLAB.

3. Place the dataset file in the same folder if available. The accepted default names are:

```text
GSE33075_expression.csv
GSE33075_expression.xlsx
gene_expression.csv
gene_expression.xlsx
expression_data.csv
expression_data.xlsx
```

4. Press **Run**.

5. Check the generated folder:

```text
Project_Output/
```

If no dataset is available, the code will still run using a demo expression dataset.

---

## Repository Structure

```text
├── MATLAB_Code/
│   └── CML_Numerical_Project.m       # Complete MATLAB workflow
│
├── Projects/
│   └── CML-Numerical-Modeling-MATLAB.md
│
├── Project_Output/                   # Created after running the MATLAB script
│   ├── Figures/
│   ├── Differential_Expression_Results.xlsx
│   ├── Treatment_Scenario_Summary.xlsx
│   └── Project_Summary.txt
│
└── README.md                         # Main portfolio page
```

---

## Key Results

The simulation and analysis demonstrate that:

- Untreated CML conditions show increased leukemic activity.
- Treatment reduces leukemic cell burden in the numerical model.
- Stronger treatment values lead to lower final leukemic cell levels.
- The mRNA expression signal decreases indirectly as leukemic burden decreases.
- Differential expression analysis can identify genes with altered expression between disease and reference groups.
- Numerical modeling and transcriptomic analysis can be combined to represent therapy response in a cancer-related biological system.

---

## Skills Demonstrated

- MATLAB programming
- Numerical analysis
- Euler’s Method
- RK4 Method
- ode45 ODE solving
- Ordinary differential equation modeling
- Gene expression data analysis
- Differential expression analysis
- Benjamini-Hochberg correction
- Scientific visualization
- Cancer modeling
- Chronic Myeloid Leukemia biology
- Computational biology documentation

---

## Professional Summary

This project demonstrates how MATLAB-based numerical methods can be used to model cancer-related molecular dynamics and targeted therapy response. By combining ODE simulation, gene expression analysis, and scientific visualization, the workflow provides a reproducible computational biology example focused on BCR-ABL1-driven Chronic Myeloid Leukemia.

---

## Author

**Shaimaa Mohamed El Haddad**  
Biomedical Science Student  
Computational Biology & Genomics Concentration

---

## Short Repository Description

MATLAB-based computational biology project modeling BCR-ABL1-related CML dynamics and targeted therapy response using Euler’s Method, RK4, ode45, and gene expression analysis.
