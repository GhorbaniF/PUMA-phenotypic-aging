PUMA: Phenotypic Unsupervised Model of Aging

PUMA is an unsupervised framework for identifying multidimensional phenotypic aging dimensions from large-scale population cohorts and evaluating their associations with future health outcomes.

The framework was developed using the Dutch Lifelines cohort and is designed to characterize heterogeneity in aging by integrating behavioral, psychological, social, physical, environmental, and biomedical phenotypes.

Features
-Preprocessing of large-scale phenotypic data.
-Construction of multidimensional phenotypic aging dimensions.
-Identification of phenotypic variables contributing to each aging dimension.
-Association analyses between phenotypic dimensions and age-related diseases.








![image](https://github.com/user-attachments/assets/bb8364cb-bcad-40e2-9593-af76cf24ef0
Workflow outlining the steps of the analysis used in this study:(I) Selecting baseline participants and variables. (II) Identifying aging-associated principal components (PCs). (III) Identifying the key baseline variables contributing to each aging-associated PC. (IV) Creating a single score from the key baseline variables contributing to each aging-associated PC referred to as “PC-composite score (PC-CS)”. (V) Evaluating the onset of age-related diseases in baseline participants. (VI) Assessing whether a defined range of PC-CSs show a statistically significant association with disease onset. (VII) Identifying critical risk factors and risk of disease onset in high-risk individuals. Image generated using BioRender.com. 


# Running the Aging Risk Factor Pipeline
This pipeline analyzes aging-related risk factors using R. It is designed to be executed on an HPC (High-Performance Computing) cluster environment.

***NOTE:*** This pipeline is optimized for execution on an HPC cluster. Make sure to load the appropriate module before running the script. The pipeline was developed and tested using R version 4.2.1, loaded via the following module on the HPC:

```
module load R/4.2.1-foss-2022a-bare
```

Once the module is loaded, run the script using:

```
Rscript aging_risk_pipeline.r
```
