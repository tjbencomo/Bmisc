context("naming")

test_that("Excel mangled gene symbols are corrected", {
  expect_equal(c("MARCH2"), check.genes(c("2-Mar"), c("chr19")))
  expect_equal(c("MARCH2", "MARC2", "2-Mar"),
               check.genes(rep("2-Mar", 3), c("chr19", "chr1", "chr3")))
})

