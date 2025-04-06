This repository contains the code files used to implement the analyses described in our paper, "**Profiling of BAL Exosomes Reveals Distinct Molecular Endotypes in Acute Respiratory Distress Syndrome.**" Below is a summary of the functionality provided by each file:

**1. **data_preprocessing_batch_effect_correction_impute.R:****
Performs data preprocessing, normalization, missing value imputation, and examines and removes batch effects.

**2. **DE.R:****
Identifies differentially expressed proteins (DEPs) between groups.

**3. **correlation.R:****
Conducts correlation analysis between DEPs and key clinical variables.

**4. **combine_correlation.R:****
Provides helper functions to combine and summarize DEPs that show significant correlations with clinical variables.

**5. clustering_correlation.R:**
Applies unsupervised K-means clustering to identify potential endotypes in ARDS patients based on correlated protein sets.

**6. endotype.R:**
Determines clinical variables that differ significantly among the identified endotypes.

**7. clustering_DE_annotation.R:**
Annotates and functionally enriches the upregulated proteins for each endotype.

**8. confounder.R:**
Examines whether age acts as a potential confounder between ARDS patients and healthy controls.

**9. figure.R:**
Contains functions to generate all the figures presented in the paper.

Each file has been developed to ensure reproducibility and clarity in the analytical workflow, facilitating further research and collaboration in the field.
