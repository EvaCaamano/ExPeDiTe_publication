# ExPeDiTe_publication

Expedite is a project with Dharani Hapangama's group at the Liverpool Women's Hospital. The main post-doc working on it is Dr Christopher Hill. Expedite has collated a number of samples of women on early stages of pregnancy that went to the emergency unit. These women provided samples of serum and urine and the team has undertaken NMR metabolomics. These women  went to have a number of different outcomes: ectopic pregnancy, pregnancy of unknown location, miscarriage or normal pregnancy.

The original aim of the study was to find markers for ectopic but prelim analysis do not see this feasible so we are aiming to find markers between normal and the rest. If the marker is validated the team could try to implement it in clinic as a way to screen women and classify them on more at risk or less at risk and divert resources to their treatment/follow up accordingly.

## Analysis caveats
Experiments were run in batches that were unbalanced. Data has been transferred to us log2 transformed and after passing through Combat. NMR signals effect sizes are small.

## Analysis plan
This is a biomarker discovery project so the strategy involved the filtering of strong signals using univariate models to pre-select variables of interest and then follow up with an attempt to find a composite marker.

Briefly, metabolite data was split in 80/20% train/validation sets. The train set was used to undergo variable selection firstly by a univariate statistical approach. The package Limma was used to generate linear mixed models to find metabolite signals different between the groups of study considering relevant metadata (BMI, Age and Gestation – these covariates were chosen as they explain a significant proportion of variance of the data) with significance adjusted for false discovery rate. Significant signals at 5% FDR were taken forward for multivariate selection. The train dataset was subjected to a 10-fold cross-validation in a 90/10 split, where each 90% split underwent 100 rounds of Least absolute shrinkage and selection operator (LASSO,glmnet package ) selection. Signals that were selected in at least 80% of the rounds and more than 8 folds were taken forward for modelling. Random forest  models were constructed with selected metabolites and metadata variables and performances were assessed in both train and validation datasets. Further tests were done by building generalised linear models and calculating AUROCs. A depiction of the workflow is shown in the image below.

![alt text](https://github.com/EvaCaamano/ExPeDiTe_publication/blob/main/WorkflowExpedite.png)

We were able to create models that predict with high accuracy  LNSP vs Other. 


## Further opportunities
As a follow up it would be relevant to study the LSNP cohort to see if there is variation within the group that can be associated to a known factor or may be indicator of something else. Some of these women went to develop pre-eclampsia or have an spontaneous pre-term birth (data not available to the CBF) and it would be worthwhile exploring potential patient subgroups. 



