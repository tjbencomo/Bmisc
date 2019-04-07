context("naming")

test_that("Excel mangled gene symbols are corrected", {
  expect_equal(c("MARCH2"), check_genes(c("2-Mar"), c("chr19")))
  expect_equal(c("MARCH2", "MARC2", "2-Mar"),
               check_genes(rep("2-Mar", 3), c("chr19", "chr1", "chr3")))
})

test_that("Correct names stay the same", {
  expect_equal(c("NONE,NONE"), check_genes(c("NONE,NONE"), c("chrM")))
  expect_equal(c("B3GALT6"), check_genes(c("B3GALT6"), c("chr1")))
})

test_that("Mitochondrial genes work", {
  expect_equal("MT-ND1", check_genes("ND1", "chrM"))
  expect_equal("IVNS1ABP", check_genes("ND1", "chr1"))
  expect_equal("ND1", check_genes("ND1", NA))
})

test_that("No name is recommended if no chromosome info", {
  expect_equal("2-Mar", check_genes("2-Mar", NA))
  expect_equal("TP53", check_genes("TP53", NA))
})

test_that("Gibberish symbols aren't corrected", {
  expect_equal("asdf", check_genes("asdf", NA))
  expect_equal("asdf", check_genes("asdf", "chr1"))

})
