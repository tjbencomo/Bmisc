#' Correct old or Excel damaged HGNC symbols
#'
#' @param original.symbol gene symbols to check
#' @param chromosome chromosomes for each gene symbol
#'
#' @return vector of predicted gene names based on chromosome info
#' @export
#'
check_genes <- function(original.symbol, chromosome) {
  mapping <- hgnc_symbol_info
  chromosome.number <- stringr::str_remove(chromosome, "chr")
  chromosome.number[chromosome.number == "M"] <- "mitochondria"
  gene.report <- HGNChelper::checkGeneSymbols(original.symbol)
  gene.report$original.chromosome <- chromosome.number

  results <- gene.report %>%
    dplyr::distinct(x, Suggested.Symbol, original.chromosome) %>%
    tidyr::separate_rows(Suggested.Symbol, sep = " /// ") %>%
    dplyr::left_join(mapping,
                     by = c("Suggested.Symbol" = "Approved symbol")) %>%
    dplyr::group_by(x) %>%
    dplyr::filter(original.chromosome == Chromosome) %>%
    dplyr::mutate(original.chromosome = dplyr::case_when(
      original.chromosome == "mitochondria" ~ "chrM",
      original.chromosome != "mitochondria" ~
        stringr::str_c("chr", original.chromosome, sep = "")
    )) %>%
    dplyr::distinct(x, Suggested.Symbol, original.chromosome)

  # if multiple rows exist for the same original name with the same chromosome
  # it is likely the old gene got split up - best to keep old name, so discard
  # these rows
  results <- results %>%
    dplyr::group_by(x) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n == length(unique(original.chromosome))) %>%
    dplyr::select(-n)

  df <- data.frame(hugo_symbol=gene.report$x,
                   chromosome=chromosome,
                   Approved=gene.report$Approved)
  df$hugo_symbol <- as.character(df$hugo_symbol)
  df$chromosome <- as.character(df$chromosome)

  report <- df %>%
    dplyr::left_join(results, by = c("hugo_symbol" = "x",
                              "chromosome" = "original.chromosome")) %>%
    dplyr::mutate(predicted.symbol = dplyr::case_when(
      Approved == T ~ hugo_symbol,
      Approved == F & !is.na(Suggested.Symbol) ~ Suggested.Symbol,
      Approved == F & is.na(Suggested.Symbol) ~ hugo_symbol
    )) %>%
    dplyr::pull(predicted.symbol)
  return(report)
}

#' Extract Amino Acid Variant
#'
#' @param aachange string containing protein sequence variant
#' @param type string indicating whether to extract the reference or alternate
#' variant. Use 'ref' to indicate reference and 'alt' to indicate the alternate
#' variant.
#'
#' @return amino acid variant specified by type
#' @export
#'
extract_aa <- function(aachange, type) {
  reference.types <- c("reference", "ref")
  alt.types <- c("mutated", "mut", "alternate", "alt")
  if (type %in% reference.types) {
    index = 1
  } else if (type %in% alt.types) {
    index = 2
  } else {
    print("Bad type!")
    return(NA)
  }
  results <- stringr::str_extract_all(aachange, '([A-Z]+[a-z]*|\\*|fs|del)', simplify = T)
  if(dim(results)[2] == 0 || (dim(results)[2] < 2 && index == 2)) {
    return(rep(NA, dim(results)[1]))
  }
  results <- results[, index]
  ifelse(results == "" | results == "UNKNOWN", NA, results)
}

#' Extract Amino Acid Variant Position
#'
#' @param aachange string containing protein sequence variant
#'
#' @return amino acid variant's position
#' @export
#'
extract_position <- function(aachange) {
  results <- stringr::str_extract(aachange, '[0-9]+_[0-9]+|[0-9]+')
  ifelse(results == "", NA, results)
}

#' Split protein sequence variant column into individual columns
#' Divides column containing protein sequence variants into unique reference,
#' position, and variant columns.
#'
#' @param df dataframe with protein sequence variant column
#' @param column name of the column to split in df
#'
#' @return updated df with reference_aa, position_aa, and alternate_aa columns
#' @export
#'
split_aachange <- function(df, column = "aachange") {
  aachanges <- df[[column]]
  ref_aa <- extract_aa(aachanges, type = "ref")
  position <- extract_position(aachanges)
  alt_aa <- extract_aa(aachanges, type = "alt")
  return(tibble::add_column(df, reference_aa = ref_aa,
                            position_aa = position,
                            alternate_aa = alt_aa))
}
