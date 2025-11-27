# Reconstruction of breast cancer subtypes associated with prognosis through multi-omics integrated analysis
![pipeline](https://github.com/user-attachments/assets/d1d2debb-327e-4a30-a962-a9f0a2df99ee)
This study integrated gene expression data, microRNA expression data, and protein expression data obtained from tumor tissues into *MOFA2* for 481 breast cancer patients in the The Cancer Genome Atlas (TCGA) cohort who had not received prior treatment or had no history of disease before their breast cancer diagnosis. The goal was to present new subgroups that better explain prognostic differences than the existing PAM50-based breast cancer subtypes and to identify the key omics signatures and biological pathways driving these differences.

## Background
### 1) Breast Cancer
![background](https://github.com/user-attachments/assets/e16e945b-851c-449f-8609-a785483f0b9a)
* Breast cancer is a highly clinically and molecularly heterogeneous cancer type, and the well-established PAM50-based subtypes are widely used to predict prognosis and treatment response.
* However, multiple studies have shown that heterogeneity exists even within the same subtype, and clinical outcomes vary significan.
* Based on these limitations, we sought to integrate multiple omics data from breast cancer to propose a new subtype that better explains prognostic differences than the existing PAM50 subtype.
  
### 2) MOFA2(Multi Omics Factor Analysis)
![mofa](https://github.com/user-attachments/assets/5e9cbab6-ab26-4fd7-a64f-00e60fd3da6c)
MOFA2 decomposes various omics data as shown in the figure.
&epsi;
## Result
### 1) MOFA2 Factor & Cox PH Model  
![mofa_result](https://github.com/user-attachments/assets/98995f06-a1fe-48aa-84e7-337fb4686931)

### 2) Cluster Analysis
![cluster_anaylsis](https://github.com/user-attachments/assets/3f5ac25b-7c0b-4195-9e08-54fa75e7518a)

### 3) Survival Analysis
![survival_anaylsis](https://github.com/user-attachments/assets/082dc5e2-4291-4c5e-a9e9-2cc593e0b15b)

### 4) Multi-omics Signature Anaylsis
![MO_signature_analysis](https://github.com/user-attachments/assets/efd6801b-65ea-4e88-886b-fe33e4429cde)

### 5) GSEA Anaylsis
![GSEA_C5](https://github.com/user-attachments/assets/8548b86f-79c8-4e89-8a5e-a96b98b96ab7)
![GSEA_C2](https://github.com/user-attachments/assets/1cfa9dd9-3334-4b06-a270-bfae7879a376)

## Conclusion
* It demonstrates that multi-omics integration can compensate for prognostic heterogeneity within subtypes that are difficult to capture using existing PAM50-based subtypes.
* The latent factors provided by MOFA2 captured common variability across omics, contributing to integrated classification.
* Signatures derived via DIABLO simultaneously identified features with high inter-omics correlations, presenting interpretable molecular networks.
* This approach may provide practical insights for subtype reclassification or developing prognostic prediction markers in clinical contexts.

