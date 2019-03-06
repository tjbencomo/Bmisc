# File: tcga_utils.R
# Author: Tomas Bencomo
# Script containing functions to interface with the
# TCGA UCSC Xena data. Easily load multiple data modalities
# for one or multiple cancers. Specify genes of interest
# or load the entire data. It is expected the data files
# follow the default naming template from Xena. 

require(tidyverse)
require(data.table)

# The naming templates used by XENA for each data type
DEFAULT.RNASEQ.FILENAME <- 'TCGA.%s.sampleMap__HiSeqV2.gz'
DEFAULT.CNV.FILENAME <- 'TCGA.%s.sampleMap__Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz'
DEFAULT.CLINICAL.FILENAME <- 'TCGA.%s.sampleMap__%s_clinicalMatrix.gz'


#' Swap axes of the dataframe. Patients are now rows and other data columns
#'
#' @param df Dataframe with patient info as columns and other data as rows
#'
#' @return Dataframe where each row is a patient and each column is some type of TCGA variable
#' @export
#'
#' @examples
transpose.data <- function(df) {
  temp <- data.frame(t(df[-1]))
  colnames(df)[1] <- 'Gene Symbol'
  colnames(temp) <- df %>% pull(`Gene Symbol`)
  setDT(temp, keep.rownames = TRUE)
  colnames(temp)[1] <- "sampleID" 
  return(temp)
}

#' Load dataframe with RNASeq data
#'
#' @param cancer TCGA study (use acronym i.e. SKCM)
#' @param data.dir Directory where the data files are stored
#' @param genes Genes to include in dataframe. NULL loads all genes
#' @param file Filename of the file containing RNASeq data
#'
#' @return Dataframe where each row is one sample and each column is
#' the expression for that gene
#' @export
#'
#' @examples
load.rnaseq <- function(cancer, data.dir, genes = NULL, file = NULL) {
  if (is.null(file)) {
    filepath.rna <- file.path(data.dir, sprintf(DEFAULT.RNASEQ.FILENAME, cancer))
  } else {
    filepath.rna <- file.path(data.dir, file)
  }
  rna.data <- read_delim(filepath.rna, "\t", 
                         escape_double = FALSE, trim_ws = TRUE, 
                         col_types = cols(), progress = FALSE)
  if(!is.null(genes)) {
    rna.data <- rna.data %>% filter(sample %in% genes)
  }
  rna.data <- transpose.data(rna.data)
  colnames(rna.data)[2:ncol(rna.data)] <- sapply(colnames(rna.data)[2:ncol(rna.data)], 
                                                 function(x) paste(x,'.rna',sep=''))
  return(rna.data)
}

#' Load dataframe with Gistic2 thresholded CNV data
#'
#' @param cancer TCGA study (use acronym i.e. SKCM)
#' @param data.dir Directory where the data files are stored
#' @param genes Genes to include in dataframe. NULL loads all genes
#' @param file Filename of the file containing CNV data
#'
#' @return Dataframe where each row is one sample and each column is
#' the Gistic2 score for that gene
#' @export
#'
#' @examples
load.cnv <- function(cancer, data.dir, genes = NULL, file = NULL) {
  if (is.null(file)) {
    filepath.cnv <- file.path(data.dir, sprintf(DEFAULT.CNV.FILENAME, cancer))
  } else {
    filepath.cnv <- file.path(data.dir, file)
  }
  cnv.data <- read_delim(filepath.cnv, "\t",
                         escape_double=FALSE, trim_ws=TRUE,
                         col_types = cols(), progress=FALSE)
  if(!is.null(genes)) {
    cnv.data <- cnv.data %>% filter(`Gene Symbol` %in% genes)
  }
  cnv.data <- transpose.data(cnv.data)
  colnames(cnv.data)[2:ncol(cnv.data)] <- sapply(colnames(cnv.data)[2:ncol(cnv.data)], 
                                                 function(x) paste(x,'.cnv',sep=''))
  return(cnv.data)
}


#' Load dataframe clinical info and optional RNASeq and CNV data
#'
#' @param cancer TCGA study (use acronym i.e. SKCM)
#' @param data.dir Directory where the data files are stored
#' @param genes Genes to include in dataframe. NULL loads all genes
#' @param include.cnv Whether to load copy number data
#' @param include.rna Whether to load rnaseq data
#' @param file Filename of the file containing clinical data
#'
#' @return Dataframe where each row is a patient sample with columns containing clinical
#' info as well as optional RNASeq/CNV data. All patients have complete Overall Survival
#' data. Solid Tissue Normal samples are not included.
#' @export
#'
#' @examples
load.cohort <- function(cancer, data.dir, genes = NULL, include.cnv = TRUE, 
                        include.rna = TRUE, file = NULL) {
  if (is.null(file)) {
    filepath.clinical <- file.path(data.dir, sprintf(DEFAULT.CLINICAL.FILENAME, cancer, cancer))
  } else {
    filepath.clinical <- file.path(data.dir, file)
  }
  cohort <- read_delim(filepath.clinical, "\t", escape_double = FALSE,
                       trim_ws = TRUE, 
                       col_types = cols(), 
                       progress = FALSE)
  if(include.rna) {
    rna.data <- load.rnaseq(cancer, data.dir, genes)
    cohort <- inner_join(cohort, rna.data, by=c('sampleID'))
  }
  if (include.cnv) {
    cnv.data <- load.cnv(cancer, data.dir, genes)
    cohort <- inner_join(cohort, cnv.data, by=c('sampleID'))
  }
  cohort <- cohort %>% filter(sample_type != 'Solid Tissue Normal')
  cohort <- drop_na(cohort, `OS`, `OS.time`)
  return(cohort)
}

#' Load multiple TCGA studies at once
#'
#' @param cancers TCGA studies (use acronym i.e. SKCM)
#' @param data.dir Directory where the data files are stored
#' @param genes Genes to include in dataframe. NULL loads all genes
#'
#' @return Dataframe where rows are patients and columns include clinical, rnaseq, and
#' cnv data. The `cohort`` column contains the TCGA study each patient belongs to.
#' @export
#'
#' @examples
load.cohorts <- function(cancers, data.dir, genes) {
  df <- data.frame()
  for(cancer in cancers) {
    current.df <- load.cohort(cancer, data.dir, genes, include.cnv = T, include.rna = T)
    current.df$cohort <- rep(cancer, dim(current.df)[1])
    df <- bind_rows(df, select(current.df, sampleID, cohort, C2orf54.rna, C2orf54.cnv, RET.rna, RET.cnv))
  }
  return(df)
}
