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

test_that("Amino acid extraction works on normal variant", {
  expect_equal("A", extract_aa("p.A430T", type = "ref"))
  expect_equal("T", extract_aa("p.A430T", type = "alt"))
  expect_equal("430", extract_position("p.A430T"))
  expect_equal(c("430", "367"), extract_position(c("p.A430T", "p.K367L")))
  expect_equal(c("A", "K"), extract_aa(c("p.A430T", "p.K367L"), type = "ref"))
  expect_equal(c("T", "L"), extract_aa(c("p.A430T", "p.K367L"), type = "alt"))
})

test_that("Extraction works for frameshift variants", {
  expect_equal("del", extract_aa("p.365_365del", type = "ref"))
  expect_equal("365_365", extract_position("p.365_365del"))
  expect_equal(NA, extract_aa("p.365_365del", type = "alt"))
  expect_equal("Q", extract_aa("p.Q357fs", type = "ref"))
  expect_equal("357", extract_position("p.Q357fs"))
  expect_equal("fs", extract_aa("p.Q357fs", type = "alt"))

})

test_that("NAs returned on empty protein sequence variants", {
  expect_equal(NA, extract_aa("p.-", type = "ref"))
  expect_equal(NA, extract_aa("p.-", type = "alt"))
  expect_equal(NA, extract_aa(".", type = "ref"))
  expect_equal(NA, extract_aa(".", type = "ref"))
  expect_equal(c(NA, NA), extract_aa(c(".", "p.-"), type = "alt"))
  expect_equal(c(NA, NA), extract_position(c(".", "p.-")))
})
