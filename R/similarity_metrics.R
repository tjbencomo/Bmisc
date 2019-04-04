library(philentropy)

jaccard.coef <- function(A, B) {
  1 - jaccard(A, B, testNA=F)
}