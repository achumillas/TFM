## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  fig.width = 8,
  fig.height = 6,
  comment   = "#>"
)

## ----install, eval = FALSE----------------------------------------------------
# 
# # Install from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ROCpAIv2")
# 
# library("ROCpAIv2")

## ----include=FALSE------------------------------------------------------------

library(dplyr)
library("ROCpAUC")


## ----Set Seed-----------------------------------------------------------------

# Set a random seed for reproducibility of all bootstrap and stochastic results
set.seed(1234567)


## ----Load Sp Data-------------------------------------------------------------

# RNA-seq data for Schizosaccharomyces pombe (fission yeast) — main demo dataset
library(fission)

# Load the fission SummarizedExperiment (RNA-seq, S. pombe, wild-type vs mut7 strain)
data(fission)

# Select five genes of interest for all package demonstrations
genes_of_interest = c(
  "SPNCRNA.1080",   # non-coding RNA, used as a low-AUC control
  "SPAC186.08c",    # protein-coding gene
  "SPNCRNA.1420",   # non-coding RNA
  "SPCC70.08c",     # protein-coding gene
  "SPAC212.04c"     # protein-coding gene
)

# Build a plain data.frame version for demos that use data.frame input.
# Rows = samples, first column = binary strain label (wt vs mut7), rest = expression.
strain_vector = colData(fission)$strain
expr_matrix_t = t(assay(fission)[genes_of_interest, ])
Sp = as.data.frame(cbind(strain = strain_vector, expr_matrix_t))

# Inspect the first rows and verify column types
head(Sp)
col_classes = sapply(Sp, class)
table(col_classes)


## ----demo-pointscurve-df------------------------------------------------------

# Input type 1: plain data.frame
# Plot the empirical ROC curve for the first gene using the pre-built Sp data.frame
pointsCurve(Sp[, 1], as.numeric(Sp[, 2]),
            plot  = TRUE,
            label = colnames(Sp)[2])


## ----demo-pointscurve-se------------------------------------------------------

# Input type 2: vectors extracted directly from a SummarizedExperiment
# Extract outcome and predictor from the SE object and pass them as vectors
outcome = SummarizedExperiment::colData(fission)$strain
expr    = SummarizedExperiment::assay(fission)["SPNCRNA.1080", ]

pointsCurve(as.numeric(outcome), as.numeric(expr),
            plot  = TRUE,
            label = "SPNCRNA.1080 (from SE)")


## ----demo-pauc-fpr-se---------------------------------------------------------

# Input type 1: SummarizedExperiment, FPR domain [0, 0.25]
pAUC(fission, low.value = 0, up.value = 0.25,
     axis = "x", selection = genes_of_interest, variable = "strain")


## ----demo-pauc-tpr-se---------------------------------------------------------

# Input type 2: SummarizedExperiment, TPR domain [0.9, 1]
pAUC(fission, low.value = 0.9, up.value = 1,
     axis = "y", selection = genes_of_interest, variable = "strain")


## ----demo-pauc-df-------------------------------------------------------------

# Input type 3: plain data.frame, FPR domain
pAUC(Sp, low.value = 0, up.value = 0.25, axis = "x")


## ----demo-mcpauc-se-tibble----------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with plots)
res_MC = MCpAUC(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  plot = TRUE, plot_type = "both", plots_per_page = 5
)
res_MC


## ----demo-mcpauc-se-se--------------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
se_MC = MCpAUC(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  output_as_SE = TRUE
)
assay(se_MC)


## ----demo-mcpauc-df-----------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
MCpAUC(Sp, low.value = 0, up.value = 0.25)


## ----demo-tpauc-se-tibble-----------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with plots)
res_Tp = TpAUC(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  plot = TRUE, plot_type = "both", plots_per_page = 5
)
res_Tp


## ----demo-tpauc-se-se---------------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
se_Tp = TpAUC(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  output_as_SE = TRUE
)
assay(se_Tp)


## ----demo-tpauc-df------------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
TpAUC(Sp, low.value = 0, up.value = 0.25)


## ----demo-npauc-se-tibble-----------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with plots)
res_Np = NpAUC(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  plot = TRUE, plot_type = "both", plots_per_page = 5
)
res_Np


## ----demo-npauc-se-se---------------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
se_Np = NpAUC(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  output_as_SE = TRUE
)
assay(se_Np)


## ----demo-npauc-df------------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
NpAUC(Sp, low.value = 0.9, up.value = 1)


## ----demo-fpauc-se-tibble-----------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with plots)
res_Fp = FpAUC(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  plot = TRUE, plot_type = "both", plots_per_page = 5
)
res_Fp


## ----demo-fpauc-se-se---------------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
se_Fp = FpAUC(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  output_as_SE = TRUE
)
assay(se_Fp)


## ----demo-fpauc-df------------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
FpAUC(Sp, low.value = 0.9, up.value = 1)


## ----demo-allindices-se-------------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with plots)
allIndices(
  fission,
  fpr.low = 0, fpr.up = 0.25,
  tpr.low = 0.9, tpr.up = 1,
  selection = genes_of_interest, variable = "strain",
  plot = TRUE, plot_type = "both", plots_per_page = 5
)


## ----demo-allindices-df-------------------------------------------------------

# Input type 2: plain data.frame → Output: tibble
allIndices(
  Sp,
  fpr.low = 0, fpr.up = 0.25,
  tpr.low = 0.9, tpr.up = 1
)


## ----demo-mcpaucboot-se-tibble------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with diagnostic plots)
boot_MC = MCpAUCboot(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  r = 100, plot = TRUE, plot_type = "both"
)
boot_MC


## ----demo-mcpaucboot-se-se----------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
boot_MC_se = MCpAUCboot(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  r = 100, output_as_SE = TRUE
)
assay(boot_MC_se)


## ----demo-mcpaucboot-df-------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
MCpAUCboot(Sp, low.value = 0, up.value = 0.25, r = 100)


## ----demo-tpaucboot-se--------------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with diagnostic plots)
boot_Tp = TpAUCboot(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  r = 100, plot = TRUE, plot_type = "both"
)
boot_Tp


## ----demo-tpaucboot-df--------------------------------------------------------

# Input type 2: plain data.frame → Output: tibble
TpAUCboot(Sp, low.value = 0, up.value = 0.25, r = 100)


## ----demo-npaucboot-se--------------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with diagnostic plots)
boot_Np = NpAUCboot(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  r = 100, plot = TRUE, plot_type = "both"
)
boot_Np


## ----demo-npaucboot-df--------------------------------------------------------

# Input type 2: plain data.frame → Output: tibble
NpAUCboot(Sp, low.value = 0.9, up.value = 1, r = 100)


## ----demo-fpaucboot-se-tibble-------------------------------------------------

# Input type 1: SummarizedExperiment → Output: tibble (with diagnostic plots)
boot_Fp = FpAUCboot(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  r = 100, plot = TRUE, plot_type = "both"
)
boot_Fp


## ----demo-fpaucboot-se-se-----------------------------------------------------

# Input type 2: SummarizedExperiment → Output: SummarizedExperiment
boot_Fp_se = FpAUCboot(
  fission,
  low.value = 0.9, up.value = 1,
  selection = genes_of_interest, variable = "strain",
  r = 100, output_as_SE = TRUE
)
assay(boot_Fp_se)


## ----demo-fpaucboot-df--------------------------------------------------------

# Input type 3: plain data.frame → Output: tibble
FpAUCboot(Sp, low.value = 0.9, up.value = 1, r = 100)


## ----summary-point-estimates--------------------------------------------------

summary_tbl = res_MC %>%
  dplyr::rename(pAUC_FPR = pAUC) %>%
  dplyr::left_join(
    dplyr::select(res_Tp, Variable, TpAUC),
    by = "Variable"
  ) %>%
  dplyr::left_join(
    dplyr::select(res_Np, Variable, NpAUC, pAUC_TPR = pAUC),
    by = "Variable"
  ) %>%
  dplyr::left_join(
    dplyr::select(res_Fp, Variable, FpAUC),
    by = "Variable"
  ) %>%
  dplyr::select(Variable, pAUC_FPR, MCpAUC, TpAUC, pAUC_TPR, NpAUC, FpAUC)

summary_tbl


## ----summary-bootstrap--------------------------------------------------------

# Build transposed bootstrap summary tables: genes as columns, metrics as rows.
# This layout makes it easy to compare the same metric across genes at a glance.

# Helper: pivot one bootstrap result to a metric × gene matrix
make_boot_matrix = function(boot_tbl, index_col) {
  # Reshape: one row per metric (estimate, bias, sd, lwr, upr), one column per gene
  metrics = c(index_col, "bias", "sd", "lwr", "upr")
  mat = t(as.matrix(boot_tbl[, metrics]))
  colnames(mat) = boot_tbl$Variable
  rownames(mat) = metrics
  return(as.data.frame(round(mat, 6)))
}

# FPR-domain bootstrap summary (MCpAUC and TpAUC) — genes as columns
cat("Bootstrap Summary — MCpAUC (FPR in [0, 0.25])\n")
make_boot_matrix(boot_MC, "MCpAUC")


## ----summary-bootstrap-Tp-----------------------------------------------------

cat("Bootstrap Summary — TpAUC (FPR in [0, 0.25])\n")
make_boot_matrix(boot_Tp, "TpAUC")


## ----summary-bootstrap-Np-----------------------------------------------------

cat("Bootstrap Summary — NpAUC (TPR in [0.9, 1])\n")
make_boot_matrix(boot_Np, "NpAUC")


## ----summary-bootstrap-Fp-----------------------------------------------------

cat("Bootstrap Summary — FpAUC (TPR in [0.9, 1])\n")
make_boot_matrix(boot_Fp, "FpAUC")


## ----demo-se-assay------------------------------------------------------------

se_MC = MCpAUC(
  fission,
  low.value = 0, up.value = 0.25,
  selection = genes_of_interest, variable = "strain",
  output_as_SE = TRUE
)

# The assay matrix: rows = metrics (MCpAUC and pAUC), columns = predictor variables
assay(se_MC)


## ----demo-se-rowdata----------------------------------------------------------

# Row metadata: names of the computed metrics stored in the assay
rowData(se_MC)


## ----demo-se-coldata----------------------------------------------------------

# Column metadata: names of the predictor variables
colData(se_MC)


