# ASD
To characterize alterations in cortical thickness-based MBNs among adolescent males with ASD and to test whether these network changes spatially correspond to normative PET-derived neurotransmitter receptor/transporter maps

# Overall Analytical Pipeline
Step 1 – Network-level analysis:
Construct cortical thickness–based morphometric brain networks (MBNs), compute nodal and global graph-theoretical metrics, and specifically analyze connectivity edges linked to the orbitofrontal cortex (OFC).

Step 2 – Linear mixed-effects modeling (LMM):
Introduce LMM to evaluate the effect of ComBat harmonization (e.g., whether site/scanner-related variance is effectively removed from network metrics).

Step 3 – Molecular-level spatial mapping:
Perform spatial correlation between network metrics and neurotransmitter receptor/transporter maps using the JuSpace toolbox. This analysis is currently not publicly shared. A placeholder folder (molecular_mapping/) is provided with the note: “Code and data for this analysis will be made available upon manuscript acceptance.”

Step 4 – Cross-level multiple comparison correction:
Apply a unified correction for multiple comparisons across all analytical levels—specifically, nodal, global, edge-wise, and molecular/neurotransmitter levels—due to the hierarchical nature of the testing framework.

Step 5 – Exploratory correlational analyses:
Examine associations between ADOS scores and (i) OFC-related connectivity edges, and (ii) nodal degree centrality of the OFC.

Step 6 – Pharmacological sensitivity analysis:
Investigate whether observed network alterations spatially align with known distributions of drug targets (e.g., receptor systems relevant to ASD pharmacotherapy).

Step 7 – Reproducibility analysis:
Assess the robustness of key findings using resampling-based approaches (e.g., bootstrap or split-half validation).

# dataset
The participants were obtained from Autism Brain Imaging Data Exchange (ABIDE) datasets (http://fcon_1000.projects.nitrc.org/indi/abide/)

# AnalysisSteps.m
AnalysisSteps.m is aimed to calculate the node properties and the correlation coefficient. This code requires the GRETNA toolbox.

# across family correction
Given the multi-level nature of the analyses, cross-level correction for multiple comparisons is required

# To perform linear mixed-effects modeling (LMM), you can use the code in the LMM_code folder.

# step2_corr_trans.m is used to calculate JUSPACE .

# The scripts for reproducibility analysis are provided in the reproducibility directory.

# Exploratory_Data_Analysis_code
We explored the association between ADOS scores and (1) the connectivity edges linked to the orbitofrontal cortex (OFC), and (2) the nodal degree centrality of the OFC.



This research was supported by Key-Area Research and Development Program of Guangdong Province (2019B030335001), National Social Science Major Project (20&ZD296), Foshan Science and Technology Innovation Projects (2420001004544), Foshan University High-level Talent Start-up Funding (CGZ07508) and Foshan University Student Academic Funding (xsjj202519zrc11 & xsjj202519zrb03).

If problems are welcome to discuss the exchange, E-mail:2020010261@m.scnu.edu.cn

# Citation
If you find this work useful in your research, please consider citing:Correspondence between morphological similarity of the left lateral orbitofrontal cortex and neurotransmitter systems in adolescent males with autism

