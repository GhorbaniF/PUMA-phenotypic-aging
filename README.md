# PUMA_pipeline
As the global population experiences accelerated aging, there is a growing interest in understanding early risk factors that have an impact on healthy aging, aiming to reduce the risk of age-related disorders. In this study, we present the first study within the LIfelines cohort, exploring the interplay between a diverse range of phenotypic and environmental risk factors and the subsequent onset of age-related diseases. Our objective is to identify early risk factors of six age-related diseases, namely cancer, diabetes, COPD, heart failure, Parkinsons and stroke. Furthermore, we aim to identify individuals who are at high-risk of disease onset to enable early interventions and improve health outcomes in those individuals.



![image](https://github.com/user-attachments/assets/bb8364cb-bcad-40e2-9593-af76cf24ef00)
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
