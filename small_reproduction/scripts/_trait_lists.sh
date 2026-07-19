#!/usr/bin/env bash
# _trait_lists.sh — curated per-dataset trait lists (sourced by analysis scripts).
# These mirror the exact gs=() arrays in the original SLURM scripts.

# 15 immune GWAS traits (Soskic / Nathan / Cano_Gamez)
IMMUNE_TRAITS=(
  "PASS_CD_deLange2017" "PASS_Celiac" "PASS_IBD_deLange2017" "PASS_Lupus"
  "PASS_Multiple_sclerosis" "PASS_Primary_biliary_cirrhosis" "PASS_Rheumatoid_Arthritis"
  "PASS_Type_1_Diabetes" "PASS_UC_deLange2017" "UKB_460K.disease_AID_ALL"
  "UKB_460K.body_HEIGHTz" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP"
  "UKB_460K.disease_RESPIRATORY_ENT" "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED"
  "UKB_460K.disease_ASTHMA_DIAGNOSED"
)

# 8 brain GWAS traits (SEA_AD / Braun)
BRAIN_TRAITS=(
  "PASS_Parkinsons23andMe_Corces2020" "PASS_Alzheimers_Jansen2019"
  "PASS_ADHD_Demontis2018" "PASS_BIP_Mullins2021"
  "PASS_Intelligence_SavageJansen2018" "PASS_Schizophrenia_Pardinas2018"
  "PASS_MDD_Howard2019" "UKB_460K.mental_NEUROTICISM"
)

# 5 microglia signatures (SEA_AD / Braun) — files in gene_sets/microglia/
MICROGLIA_SIGS=( "HM_gs" "DAM_gs" "CRM_gs" "IRM_gs" "HLA_gs" )

# 52 T-cell phenotype signatures (Soskic / Nathan / Cano_Gamez) — gene_sets/t_cell_pheno/
TCELL_SIGS=(
  "Doublet-Myeloid" "Metallothionein" "Translation" "IL10-IL19" "OX40-EBI3"
  "CD172a-MERTK" "Th2-Activated" "CD4-CM" "CD8-Trm" "TIMD4-TIM3" "BCL2-FAM13A"
  "CellCycle-G2M" "IEG" "MAIT" "SOX4-TOX2" "Doublet-Fibroblast" "NME1-FABP5"
  "Th17-Activated" "IEG3" "RGCC-MYADM" "CD8-Naive" "Exhaustion" "ISG"
  "Th2-Resting" "CD4-Naive" "CellCycle-Late-S" "Treg" "Tph" "Doublet-Plasmablast"
  "TEMRA" "Cytotoxic" "Doublet-Platelet" "CD40LG-TXNIP" "CellCycle-S" "Th1-Like"
  "Mito" "Tfh-2" "gdT" "Poor-Quality" "Doublet-RBC" "HLA" "Heatshock" "IEG2"
  "Cytoskeleton" "CTLA4-CD38" "Doublet-Bcell" "Multi-Cytokine" "Th17-Resting"
  "Th22" "ICOS-CD38" "Tfh-1" "CD8-EM"
)
