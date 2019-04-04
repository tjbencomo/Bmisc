library(HGNChelper)
library(tidyverse)

# Returns vector of "correct names"
checkGenes <- function(original.symbol, chromosome) {
  mapping <- read_delim("~/Bmisc/results.txt", "\t", escape_double = FALSE, 
                        trim_ws = TRUE, progress = F)
  chromosome.number <- as.integer(str_extract(chromosome, "([0-9])+"))
  gene.report <- checkGeneSymbols(original.symbol)
  gene.report$original.chromosome <- chromosome.number
  
  results <- gene.report %>%
    distinct(x, Suggested.Symbol, original.chromosome) %>%
    separate_rows(Suggested.Symbol, sep = " /// ") %>%
    left_join(mapping, by = c("Suggested.Symbol" = "Approved symbol")) %>%
    group_by(x) %>%
    filter(original.chromosome == Chromosome) %>%
    mutate(original.chromosome = paste("chr", original.chromosome, sep='')) %>%
    distinct(x, Suggested.Symbol, original.chromosome)
  
  # if multiple rows exist for the same original name with the same chromosome
  # it is likely the old gene got split up - best to keep old name, so discard
  # these rows
  results <- results %>%
    group_by(x) %>%
    mutate(n = n()) %>%
    filter(n == length(unique(original.chromosome))) %>%
    select(-n)
  
  df <- data.frame(hugo_symbol=gene.report$x, 
                   chromosome=chromosome, 
                   Approved=gene.report$Approved)
  df$hugo_symbol <- as.character(df$hugo_symbol)
  df$chromosome <- as.character(df$chromosome)
  
  report <- df %>% 
    left_join(results, by = c("hugo_symbol" = "x", 
                              "chromosome" = "original.chromosome")) %>%
    mutate(predicted.symbol = case_when(
      Approved == T ~ hugo_symbol,
      Approved == F & !is.na(Suggested.Symbol) ~ Suggested.Symbol,
      Approved == F & is.na(Suggested.Symbol) ~ hugo_symbol
    )) %>% 
    pull(predicted.symbol)
  return(report)
}
