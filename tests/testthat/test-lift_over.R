context("Test lift over")

library(S4Vectors)

test_that("convert methCall to GRanges works", {
  methCallObj <- new("methCall", data = DataFrame(chr = Rle("chr1"), position =Rle(226061851), strand = Rle("+"), met_reads = Rle(1), coverage = Rle(1)))
  grObj <- as_GRanges(methCallObj)
  expect_equal(grObj@seqnames@values, factor("chr1"))
  expect_equal(grObj@ranges@start, 226061851)
  expect_equal(grObj@strand@values, factor("+", levels = c("+", "-", "*")))
})

test_that("Liftover from mm9 to mm10 works", {
  methCallObj <- new("methCall", data = DataFrame(chr = Rle(c("30", "3")), position =Rle(c(0, 133126641)), strand = Rle(rep("+", 2)), met_reads = Rle(c(1, 2)), coverage = Rle(rep(1, 2))))
  methCallObj_mapped <- map_coord(methCallObj, system.file("extdata/mm9ToMm10.over.chain", package = "scDNAmClock"))
  expect_equal(n_row(methCallObj_mapped), 1)
  expect_equal(as.character(methCallObj_mapped@data$chr@values), "3")
  expect_equal(methCallObj_mapped@data$position@values, 133463677)
})

