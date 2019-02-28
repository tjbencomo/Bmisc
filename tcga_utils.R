require(pacman)
pacman::p_load(tidyverse, data.table)

transpose.data <- function(df) {
  temp <- data.frame(t(df[-1]))
  colnames(df)[1] <- 'Gene Symbol'
  colnames(temp) <- df %>% pull(`Gene Symbol`)
  setDT(temp, keep.rownames = TRUE)
  colnames(temp)[1] <- "sampleID" 
  return(temp)
}

load.rnaseq <- function(cancer, data.dir, genes = NULL) {
  filepath.rna <- file.path(data.dir, 
                            paste('TCGA.', cancer, 
                                  '.sampleMap__HiSeqV2.gz',sep=''))
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

load.cnv <- function(cancer, data.dir, genes = NULL) {
  path.template <- '.sampleMap__Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz'
  filepath.cnv <- file.path(data.dir, 
                            paste('TCGA.',cancer, path.template, sep = ''))
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


load.cohort <- function(cancer, data.dir, genes = NULL, include.cnv = TRUE, include.rna = TRUE) {
  filepath.clinical <- file.path(data.dir, 
                                 paste('TCGA.', cancer, '.sampleMap__',cancer,
                                       '_clinicalMatrix.gz', sep = ''))
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

load.cohorts <- function(cancers, data.dir, genes) {
  df <- data.frame()
  for(cancer in cancers) {
    # print(cancer)
    current.df <- load.cohort(cancer, data.dir, genes, include.cnv = T, include.rna = T)
    current.df$cohort <- rep(cancer, dim(current.df)[1])
    df <- bind_rows(df, select(current.df, sampleID, cohort, C2orf54.rna, C2orf54.cnv, RET.rna, RET.cnv))
  }
  return(df)
}

