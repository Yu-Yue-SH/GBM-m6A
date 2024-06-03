# cellchat analysis for different m6A levels
# function: run_cellchat(data_input, meta, group, w)
# function: cellchat_comparison(cellchat_hi, cellchat_lo)


rm(list = ls())


library(Seurat)
library(CellChat)
library(patchwork)


options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)


# load data
GBM_core_stratified <- readRDS('outputs/seu/GBM_core_stratified_04.rds')


# m6a
writers <- c('METTL3', 'METTL5', 'METTL14', 'METTL16', 'WTAP', 'RBM15',
             'RBM15B', 'ZNF217', 'ZC3H13', 'VIRMA', 'CBLL1')
erasers <- c('FTO', 'ALKBH5')
readers <- c('YTHDF1', 'YTHDF2', 'YTHDF3', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3',
             'YTHDC1', 'YTHDC2', 'ELAVL1', 'HNRNPC', 'RBMX', 'HNRNPA2B1')
m6a <- c(writers, erasers, readers)


# function: run_cell_chat ####
run_cellchat <- function(data_input, meta, group = 'new_broad_type', w = 8) {
  # create object
  cellchat <- createCellChat(data_input, meta, group.by = group)

  # set database
  CellChatDB <- CellChatDB.human
  CellChatDB_use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB_use

  # preprocess
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = w)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # run
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}


# run run_cell_chat
for (i in m6a) {
  # first high then low
  for (j in c('High', 'Low')) {
    seu_obj <- subset(
      GBM_core_stratified,
      cells = which(
        GBM_core_stratified[['new_broad_type']] != 'Cancer cell' |
          GBM_core_stratified[[paste0(i, '_level')]] == j
      )
    )
    data_input <- GetAssayData(seu_obj, layer = 'data')
    meta <- seu_obj@meta.data
    cellchat <- run_cellchat(data_input, meta, 'new_broad_type', 8)
    saveRDS(
      cellchat,
      paste0('outputs/cellchat/04/raw/cellchat_', i, '_', j, '.rds')
    )
  }
}


# function: cellchat_comparison ####
cellchat_comparison <- function (cellchat_hi, cellchat_lo) {
  # merge object
  object_list <- list(lo = cellchat_lo, hi = cellchat_hi)

  cellchat_diff <- mergeCellChat(
    object_list,
    add.names = names(object_list),
    cell.prefix = TRUE
  )

  return(cellchat_diff)
}


# run cellchat_comparison
for (i in m6a) {
  # load data
  cellchat_hi <-
    readRDS(paste0('outputs/cellchat/04/raw/cellchat_', i, '_', 'High', '.rds'))
  cellchat_lo <-
    readRDS(paste0('outputs/cellchat/04/raw/cellchat_', i, '_', 'Low', '.rds'))

  # run cellchat_comparison
  cellchat_diff <- cellchat_comparison(cellchat_hi, cellchat_lo)

  saveRDS(
    cellchat_diff,
    paste0('outputs/cellchat/04/diff/cellchat_', i, '_diff.rds')
  )
}
