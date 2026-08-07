# Hepatoblastoma_CosMx-Visium

Reproducible computational workflow accompanying the manuscript:

**Spatially resolved tumor-microenvironment architecture defines therapeutic response in unresectable hepatoblastoma**

This repository contains the computational workflow used to analyze **10x Genomics Visium** spatial transcriptomics and **NanoString CosMx Spatial Molecular Imaging (SMI)** data. The pipeline integrates both modalities to characterize spatial tumor ecosystems, identify cellular neighborhoods, infer cell-cell communication, and validate relapse-associated transcriptional signatures.

---

# Repository Structure

| Directory | Description |
|-----------|-------------|
| `01_Visium` | Visium preprocessing, quality control, deconvolution and region annotation |
| `02_CosMx` | CosMx preprocessing, clustering and cell-state annotation |
| `03_Integration` | Integration of Visium and CosMx datasets |
| `04_Downstream_Analysis` | Spatial neighborhood, CellChat, pathway and network analyses |
| `05_Validation` | Relapse signature validation |
| `06_Manuscript_Figures` | Scripts for reproducing manuscript figures |

---

# Analysis Workflow

```text
Visium
      │
      ▼
Quality Control
      │
      ▼
STdeconvolve
      │
      ▼
Spatial Region Annotation
      │
      ▼
Visium–CosMx Integration
      │
      ├── CellChat
      ├── Giotto
      ├── Pathway Analysis
      ├── PPI / TF Networks
      └── Signature Validation
```

---

# Data

Visium and CosMx datasets have been deposited to GEO.

**Accession numbers:** Pending publication.

---

# Contact

**Zahra Mesrizadeh, PhD**

Department of Bioengineering  
University of California San Diego

Email: zmesriza@ucsd.edu
