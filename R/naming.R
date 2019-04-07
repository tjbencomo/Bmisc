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
  print(chromosome.number)
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
