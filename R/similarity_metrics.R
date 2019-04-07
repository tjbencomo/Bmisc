library(philentropy)

jaccard_coef <- function(A, B) {
  1 - jaccard(A, B, testNA=F)
}
