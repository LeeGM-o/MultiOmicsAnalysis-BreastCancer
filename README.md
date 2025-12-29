# Reconstruction of breast cancer subtypes associated with prognosis through multi-omics integrated analysis
![pipeline](https://github.com/user-attachments/assets/d1d2debb-327e-4a30-a962-a9f0a2df99ee)
This study integrated gene expression data, microRNA expression data, and protein expression data obtained from tumor tissues into *MOFA2* for 481 breast cancer patients in the The Cancer Genome Atlas (TCGA) cohort who had not received prior treatment or had no history of disease before their breast cancer diagnosis. The goal was to present new subgroups that better explain prognostic differences than the existing PAM50-based breast cancer subtypes and to identify the key omics signatures and biological pathways driving these differences.

## Background
### 1) Breast Cancer
![background](https://github.com/user-attachments/assets/e16e945b-851c-449f-8609-a785483f0b9a)
Breast cancer is a highly clinically and molecularly heterogeneous cancer type, and the well-established PAM50-based subtypes are widely used to predict prognosis and treatment response. However, multiple studies have shown that heterogeneity exists even within the same subtype, and clinical outcomes vary significan. Based on these limitations, we sought to integrate multiple omics data from breast cancer to propose a new subtype that better explains prognostic differences than the existing PAM50 subtype.
  
### 2) MOFA2(Multi Omics Factor Analysis)
MOFA2(Multi Omics Factor Analysis) is an unsupervised learning data integration framework that captures the variability across multiple omics data from the same sample to identify *latent factors*. These latent factors serve as compressed, low-dimensional interpretable representations of the multi-omics dataset, enabling subsequent analyses such as *clustering* and *feature selection*. Each factor assigns weights to samples and features, representing complex data as a few large, easily understandable axes.
  
![mofa](https://github.com/user-attachments/assets/5e9cbab6-ab26-4fd7-a64f-00e60fd3da6c)  

MOFA2 decomposes various omics data as shown in the figure. Since *Z* and *W* can't accurately reproduce the actual observed values, ε is added to model them as probabilistic noise. MOFA2 performs inference through **stochastic variational Bayesian** methods, which assume prior joint distributions for unobserved variables within the model.
* Citation: [Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, et al. Multi-Omics Factor Analysis-a framework for unsupervised integration of multi-omics data sets. Mol Syst Biol. 2018;14:e8124](https://www.embopress.org/doi/full/10.15252/msb.20178124)
  
## Result
### 1) MOFA2 Factor & Cox PH Model  
![mofa_result](https://github.com/user-attachments/assets/98995f06-a1fe-48aa-84e7-337fb4686931)
* MOFA2 captures the patterns and variability of each omics dataset, compressing them into 15 latent factors and presenting the results as a heatmap and bar plot.
* The heatmap shows the proportion of *Z* * *W* variance explained by each factor across all datasets, while the bar plot displays the proportion of total variance explained by each factor for each omics type.
* RNA-seq was confirmed to explain the largest variance. Furthermore, the three omics datasets collectively explain a significant portion of variance within each factor, indicating the presence of **shared biological signals across the multi-omics data**.
* To assess the prognostic value of MOFA2, a multivariate Cox proportional hazards model was applied using 15 factors. The forest plot showing hazard ratios and confidence intervals revealed that 87 events occurred among the 481 patients, confirming that **factor 3, factor 6, and factor 13** were significantly associated with survival(p<0.05).

### 2) Cluster Analysis
![cluster_anaylsis](https://github.com/user-attachments/assets/3f5ac25b-7c0b-4195-9e08-54fa75e7518a)
* To create new clusters using information from factors 3, 6, and 13, which were significantly associated with survival, we used the gap statistic and silhouette score to select the optimal k before performing cluster analysis. We selected k=3, which yielded the highest gap statistic, and k=3 also showed a high silhouette score.
* The three clusters were designated as **MOG (Multi Omics Group) 1, 2, and 3**. Dimension reduction was performed using umap to visualize the clusters in two dimensions and confirm their patterns.
* First, the existing PAM50 subtypes showed unclear clustering patterns, while the groups divided into MOG1, 2, and 3 formed distinct clusters to some extent. The classification into MOG1 and MOG2 was particularly clear.
* The bar graph confirmed that PAM50 subtypes were distributed across MOG1, MOG2, and MOG3. However, the Basal subtype had the highest proportion in MOG2, while the Luminal A subtype constituted the majority in MOG1.

### 3) Survival Analysis
![survival_anaylsis](https://github.com/user-attachments/assets/082dc5e2-4291-4c5e-a9e9-2cc593e0b15b)
* We then performed survival analysis to determine if there were differences in prognosis among the MOG groups. For comparison, we also conducted survival analysis for the PAM50 subtype.
* Among the three groups, MOG1 showed the best prognosis.
* When performing pairwise statistical tests, we found significant prognostic differences between the MOG1-2 and MOG1-3 pairs, whereas the existing PAM50 subtype showed weaker prognostic association.
* **These results demonstrate that multi-omics integration can complement and explain prognostic differences that are difficult to account for using only the existing PAM50-based subtypes, emphasizing the necessity of multi-omics integration**. 
* Subsequent analyses focused on MOG1-2 and MOG1-3, where significant prognostic differences were evident in survival curves.

### 4) Multi-omics Signature Analysis
#### 4-1. mixOmics - DIABLO
![DIABLO](https://github.com/user-attachments/assets/883fa656-4e14-420e-bdd3-9c61432e2e01)
* We employed the *DIABLO* framework from mixOmics to identify key features exhibiting expression differences between the two groups showing differential prognosis.
* *DIABLO* is a supervised learning-based framework that learns covariate relationships across different omics datasets to identify features that effectively distinguish specific groups.
> 1. When the input X is provided as a list of each omics data frame and the group class labels are provided in Y, training proceeds using the Sparse Partial Least Squares Discriminant Analysis(sPLS-DA) algorithm. This algorithm creates a weighted sum of the features from each omics in X to maximize the covariance with Y.
> 2. This process defines common latent components, similar to dimension reduction, selecting key features most contributing to group differentiation across each omics while also accounting for correlations between omics.
> 3. During training, 10-fold cross-validation was used to determine the optimal number of features for each omics dataset.
#### 4-2. Circular & Network plot
![MOS_analysis](https://github.com/user-attachments/assets/58a1ed20-ee42-4da9-b4cf-ccd9128c9704)
* The circular plot allows for an intuitive visualization of which features show high expression in specific groups and their associations with other omics features.
* In a network plot represented solely by connections between features with Pearson correlation coefficients of |0.7| or higher, a network structure showing positive correlations formed between *LCK*(protein), *hsa-miR-150*(miRNA), and *Immunoglobulin(IG) family genes*(RNA-seq), with high expression of these features observed in MOG1.
* These features were indeed reported in the paper to be associated with immune responses in the **tumor microenvironment**.   
> **LCK**: A key switch in T cell immunity, serving as an indicator reflecting the activation and TCR signaling levels of tumor-infiltrating T cells.  
> **Immunoglobulin family genes**: Reflect the degree of B cell and plasma cell infiltration.  
> **hsa-mir-150**: An immune-related miRNA specifically expressed in T cells and B cells, delivered to cells via exosomes, acting as a regulator associated with lymphocyte differentiation and maintenance of immune homeostasis.
* **These results demonstrate the characteristics of multi-omics integration, which is difficult to capture with a single omics approach. They suggest that an immune-active tumor microenvironment, where T cells and B cells are simultaneously activated in MOG1, may be associated with a favorable prognosis due to a high immune response. Furthermore, the connections between features revealed in the network structure indicate that they could serve as clinically useful indicators in multi-omics integration studies related to breast cancer prognosis.**

  
### 5) GSEA Analysis
#### 5-1. C5 - Gene Ontology Biological Process
![GSEA_C5](https://github.com/user-attachments/assets/8548b86f-79c8-4e89-8a5e-a96b98b96ab7)
#### 5-2. C2 - KEGG Pathway
![GSEA_C2](https://github.com/user-attachments/assets/1cfa9dd9-3334-4b06-a270-bfae7879a376)
* To identify which pathways were active in the two groups showing different prognoses, GSEA analysis was performed using the C5 and C2 databases from MsigDB. Enrichment analysis revealed that cell proliferation-related pathways were enriched in MOG2 and MOG3, while immune-related pathways were enriched in MOG1.
* The enrichment of immune-related pathways in MOG1 functionally supports that the network identified in the prior multi-omics signature analysis indeed reflects an immune-active tumor microenvironment.
* The enrichment of cell cycle and mitosis pathways in MOG2 and MOG3 functionally confirms that the absence of common networks observed in DIABLO and the difference in PAM50 ratios among samples are driven by a proliferation-centric tumor intrinsic program. These biological characteristics align with existing reports associating low immune infiltration with poor prognosis and provide a mechanistic basis explaining the poor clinical outcomes observed for MOG2 and MOG3 in this study.

## Conclusion
* It demonstrates that multi-omics integration can compensate for prognostic heterogeneity within subtypes that are difficult to capture using existing PAM50-based subtypes.
* The latent factors provided by MOFA2 captured common variability across omics, contributing to integrated classification.
* Signatures derived via DIABLO simultaneously identified features with high inter-omics correlations, presenting interpretable molecular networks.
* This approach may provide practical insights for subtype reclassification or developing prognostic prediction markers in clinical contexts.

## Data Availability
* The data used in this study can be downloaded from [UCSC Xena](https://xenabrowser.net/datapages/).
* Additionally, PAM50 subtype information can be obtained from the R package *TCGABiolinks*.
