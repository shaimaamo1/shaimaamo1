# MATLAB-Based Numerical Simulation of BCR-ABL1 Expression Dynamics and Targeted Therapy Response in Chronic Myeloid Leukemia

## Project Overview

This project presents a MATLAB-based computational and numerical simulation of **BCR-ABL1 expression dynamics** in **Chronic Myeloid Leukemia (CML)**. The project models the behavior of BCR-ABL1 mRNA and protein/activity levels under three biological conditions: normal cells, untreated CML cells, and treated CML cells.

The biological focus of the project is based on the role of the **BCR-ABL1 gene** in CML, where its abnormal tyrosine kinase activity promotes excessive survival and growth of myeloid cells. The project also incorporates a treatment-response scenario using **tyrosine kinase inhibitor-like inhibition**, inspired by imatinib treatment data from the public GEO dataset **GSE33075**.

## Project Aim

The main aim of this project was to model **BCR-ABL1 mRNA and protein/activity dynamics** in:

- Normal cells
- Untreated CML cells
- Treated CML cells

The expected biological behavior was that normal cells would show the lowest BCR-ABL1 activity, untreated CML cells would show the highest activity, and treated CML cells would show reduced activity after therapy.

## Mathematical Model

The model uses two main biological variables:

- **M(t):** BCR-ABL1 mRNA level
- **P(t):** BCR-ABL1 protein/activity level

The model assumes that mRNA dynamics depend on synthesis and decay, while protein/activity dynamics depend on synthesis, decay, and treatment inhibition. The treatment effect was represented using an inhibition term affecting the protein/activity level.

## Numerical Methods

The project applied two numerical simulation approaches in MATLAB:

- **Euler’s Method** for step-by-step approximation of the model behavior
- **MATLAB ode45** for solving the ordinary differential equation system

The simulation was performed over time to compare BCR-ABL1 mRNA and protein/activity levels across the three biological scenarios.

## Transcriptomic Data Analysis

The project also included analysis of the public GEO transcriptomic dataset **GSE33075**. The dataset contains samples from healthy controls, CML patients before treatment, and CML patients after one month of imatinib treatment. Expression values were imported into MATLAB, grouped by condition, and compared to identify genes that decreased after treatment.

## Key Results

The simulation showed that:

- Normal cells had low BCR-ABL1 mRNA and protein/activity levels.
- Untreated CML cells showed the highest BCR-ABL1 protein/activity levels.
- Treated CML cells showed reduced protein/activity due to the treatment inhibition term.
- Increasing treatment strength led to progressively lower BCR-ABL1 activity.
- GEO transcriptomic analysis identified downregulated genes after treatment, supporting the treatment-response scenario.

## Tools and Skills Used

- MATLAB
- Numerical Analysis
- Euler’s Method
- ode45 ODE Solver
- Gene Expression Analysis
- GEO Dataset Analysis
- Mathematical Modeling
- Biomedical Simulation
- Chronic Myeloid Leukemia Biology
- BCR-ABL1 Expression Dynamics
- Treatment Response Modeling
- Scientific Data Visualization

## Conclusion

This project demonstrated how numerical methods and computational biology can be integrated to study cancer-related gene expression behavior. By combining mathematical modeling, MATLAB simulation, and transcriptomic data analysis, the project provided a simplified computational representation of BCR-ABL1 dynamics in CML and showed how targeted therapy may reduce leukemic activity over time.

## Repository Structure

```text
├── MATLAB_Code/          # MATLAB scripts for model simulation and GEO analysis
├── Figures/              # Generated plots and simulation outputs
├── Dataset/              # Processed GEO dataset files
├── Presentation/         # Final project presentation
├── Report/               # Project report/documentation
└── README.md             # Project overview and documentation
```

## Short Repository Description

MATLAB-based numerical simulation of BCR-ABL1 mRNA and protein/activity dynamics in Chronic Myeloid Leukemia, integrating Euler’s Method, ode45, and GEO gene expression analysis.

## Suggested Repository Name

```text
CML-BCRABL1-MATLAB-Simulation
```

or

```text
BCRABL1-CML-Numerical-Modeling
```

## Author

**Shaimaa Mohamed El Haddad**  
Biomedical Science Student  
Computational Biology & Genomics Concentration

---

## Professional GitHub Description

A computational biology project modeling BCR-ABL1 dynamics and treatment response in CML using MATLAB numerical methods and GEO transcriptomic data.
