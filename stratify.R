# preprocess and stratify m6A levels in GBM_core dataset
# function: stratify_m6a_levels_in_platforms(obj, platforms, thresholds, m6a)


rm(list = ls())


library(Seurat)


# load data
GBM_core <- readRDS('data/GBmap/GBM_core.rds')


# preprocess
# change gene ids into names
counts <- GetAssayData(GBM_core, layer = 'counts')
data <- GetAssayData(GBM_core, layer = 'data')
meta_features <- GBM_core[['RNA']]@meta.features
meta_features[['ENSG_id']] <- rownames(meta_features)
rownames(meta_features) <- meta_features[['feature_name']]
rownames(counts) <- meta_features[['feature_name']]
rownames(data) <- meta_features[['feature_name']]
meta_data <- GBM_core@meta.data
re_umap <- GBM_core@reductions$umap

GBM_core <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_data,
  project = 'GBM_core'
)

GBM_core <- SetAssayData(GBM_core, layer = 'data', data)
GBM_core[['RNA']]@meta.data <- meta_features
GBM_core@reductions$umap <- re_umap


# change platform ids into names
GBM_core[['platform']] <- NA
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0008931'] <-
  'Smart-seq2'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0008953'] <-
  'STRT-seq'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0009899'] <-
  '10x 3 v2'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0009922'] <-
  '10x 3 v3'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0010010'] <-
  'CEL-seq2'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0011025'] <-
  '10x 5 v1'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0030002'] <-
  'Microwell-seq'
GBM_core[['platform']][GBM_core[['assay_ontology_term_id']] == 'EFO:0700003'] <-
  'BD Rhapsody'

saveRDS(GBM_core, 'outputs/seu/GBM_core_pp.rds')
GBM_core <- readRDS('outputs/seu/GBM_core_pp.rds')


# m6a
# 1	Jiang, X. et al. The role of m6A modification in the biological functions
# and diseases. Signal Transduct Target Ther 6, 74 (2021).
# https://doi.org/10.1038/s41392-020-00450-x
# *2	Fang, Z. et al. Role of m6A writers, erasers and readers in cancer. Exp
# Hematol Oncol 11, 45 (2022). https://doi.org/10.1186/s40164-022-00298-7
# 3	Zaccara, S., Ries, R. J. & Jaffrey, S. R. Reading, writing and erasing mRNA
# methylation. Nat Rev Mol Cell Biol 20, 608-624 (2019).
# https://doi.org/10.1038/s41580-019-0168-5
writers <- c('METTL3', 'METTL5', 'METTL14', 'METTL16', 'WTAP', 'RBM15',
             'RBM15B', 'ZNF217', 'ZC3H13', 'VIRMA', 'CBLL1')
#  VIRMA == KIAA1429, CBLL1 == HAKAI
erasers <- c('FTO', 'ALKBH5')
readers <- c('YTHDF1', 'YTHDF2', 'YTHDF3', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3',
             'YTHDC1', 'YTHDC2', 'ELAVL1', 'HNRNPC', 'RBMX', 'HNRNPA2B1')
# RBMX==HNRNPG
# all HNRNPs, while only important genes are selected
# rownames(GBM_core)[grep('^HNRNP', rownames(GBM_core))]
# EIF3s are not included
m6a <- c(writers, erasers, readers)


# function: stratify_m6a_levels_in_platforms ####
stratify_m6a_levels_in_platforms <- function (
  obj,
  platforms = c('10x 3 v2', '10x 3 v3', '10x 5 v1', 'BD Rhapsody', 'CEL-seq2',
                'Microwell-seq', 'Smart-seq2', 'STRT-seq'),
  thresholds = c(0, 0.4, 0.6, 1),
  m6a = c('METTL3', 'METTL5', 'METTL14', 'METTL16', 'WTAP', 'RBM15',
          'RBM15B', 'ZNF217', 'ZC3H13', 'VIRMA', 'CBLL1', 'FTO', 'ALKBH5',
          'YTHDF1', 'YTHDF2', 'YTHDF3', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3',
          'YTHDC1', 'YTHDC2', 'ELAVL1', 'HNRNPC', 'RBMX', 'HNRNPA2B1')
) {
  # extract data
  cancer_cells <- subset(obj, annotation_level_1 == 'Neoplastic')

  data <- as.data.frame(GetAssayData(cancer_cells, layer = 'data'))
  data <- data[m6a, ]
  data <- t(data)
  data <- cbind(data, cancer_cells@meta.data)

  # stratify m6a levels
  for (i in m6a) {
    data[paste0(i, '_level')] <- NA
    data[paste0(i, '_level')][data[[i]] == 0, ] <- 'Not detected'
    for (j in platforms) {
      cells <- data[[i]] != 0 & data[['platform']] == j
      if (sum(cells) == 0) {next}
      data[cells, paste0(i, '_level')] <- as.character(
        cut(
          data[cells, i],
          breaks = quantile(data[cells, i], thresholds),
          labels = c('Low', 'Medium', 'High'),
          include.lowest = TRUE
        )
      )
    }
  }

  # add to seurat object
  obj[['new_broad_type']] <- obj[['annotation_level_2']]
  obj@meta.data[['new_broad_type']] <-
    as.character(obj@meta.data[['new_broad_type']])
  obj[['new_broad_type']][obj[['annotation_level_1']] == 'Neoplastic'] <-
    'Cancer cell'

  for (i in m6a) {
    cancer_cells[[paste0(i, '_level')]] <- NA
    cancer_cells[[paste0(i, '_level')]]<- data[[paste0(i, '_level')]]
    obj[[paste0(i, '_level')]] <- obj[['new_broad_type']]
    obj[[paste0(i, '_level')]] <- cancer_cells[[paste0(i, '_level')]]
  }

  return(obj)
}


# run stratify_m6a_levels_in_platforms
# due to insufficient sequencing depth, CEL-seq2 is removed
GBM_core_stratified <- subset(GBM_core, platform != 'CEL-seq2')
GBM_core_stratified <- stratify_m6a_levels_in_platforms(
  GBM_core_stratified,
  platforms = c('10x 3 v2', '10x 3 v3', '10x 5 v1', 'BD Rhapsody',
                'Microwell-seq', 'Smart-seq2', 'STRT-seq'),
  thresholds = c(0, 0.4, 0.6, 1)
)

saveRDS(GBM_core_stratified, 'outputs/seu/GBM_core_stratified_04.rds')
GBM_core_stratified <- readRDS('outputs/seu/GBM_core_stratified_04.rds')
