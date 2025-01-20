# MNMI_AD

This GitHub repository holds the data and code that support the findings of the MNMI-AD study. The detailed method and major results can be found in the following report:

     * Li G, Hsu LM, Wu Y, Bozoki AC, Shih YI, Yap PT (2025) Revealing excitation-inhibition imbalance in Alzheimer's disease using multiscale neural model inversion of resting-state functional MRI. Commun Med (Lond). 5(1):17. doi: 10.1038/s43856-025-00736-7. PMID: 39814858; PMCID: PMC11735810.

The “SUBJECT” folder stores the demographic and clinical information of participants.

The “DATA” folder stores the functional connectivity (FC) and BOLD time series for 48 normal control (NC), 48 mild cognitive impairment (MCI), and 48 Alzheimer’s disease (AD) subjects. It also stores the structural connectivity (SC) for 100 HCP subjects.

The “MNMI” folder stores the Matlab codes to estimate connection parameters of individual subjects. The code was run in parallel with 24 cores in a high-performance UNC Linux computing cluster operating with Red Hat Enterprise Linux 7.

The “OUTPUT” folder stores the output (i.e., estimated connection parameters) of individual subjects.

The “ANALYSIS” folder stores the Matlab code to analyze the output of the MNMI model (i.e., statistical comparison among NC, MCI, and AD groups).

For questions, please contact Guoshi Li (guoshi_li@med.unc.edu). 
