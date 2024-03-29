# We should be able to reproduce identical plots by setting the random seed.
#
#   ggplot(...) + geom_text_repel(..., seed = 1)
#
context("seed")

library(grid)

pos_df <- function(pos) {
  data.frame(
    x = sapply(pos, function(x) {
      convertWidth(x[["x"]], "native")
    }),
    y = sapply(pos, function(x) {
      convertHeight(x[["y"]], "native")
    })
  )
}

test_that("calling geom_text_repel without seed creates different plots", {

  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # Make a plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
  png(png_file)
  p1 <- ggplot(dat1) + geom_text_repel(aes(wt, mpg, label = label))
  print(p1)
  grid.force()
  pos1 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()
  unlink("testthat_test-seed1.png")

  # Make a second plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed2")
  png(png_file)
  p2 <- ggplot(dat1) + geom_text_repel(aes(wt, mpg, label = label))
  print(p2)
  grid.force()
  pos2 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Confirm that the label positions are identical.
  expect_true(nrow(pos1) == nrow(dat1))
  expect_true(nrow(pos2) == nrow(dat1))
  expect_true(!identical(pos1, pos2))
})

test_that("calling geom_text_repel with seed creates identical plots", {
  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # Make a plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
  png(png_file)
  p1 <- ggplot(dat1) + geom_text_repel(aes(wt, mpg, label = label), seed = 10)
  print(p1)
  grid.force()
  pos1 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Make a second plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed2")
  png(png_file)
  p2 <- ggplot(dat1) + geom_text_repel(aes(wt, mpg, label = label), seed = 10)
  print(p2)
  grid.force()
  pos2 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Confirm that the label positions are identical.
  expect_true(nrow(pos1) == nrow(dat1))
  expect_true(nrow(pos2) == nrow(dat1))
  expect_equal(pos1, pos2, tolerance = 1e-3)
})

test_that("calling geom_text_repel without seed does not remove entropy", {
  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # One random number will be generated after each plot and accumulated
  random_seq = c()
  for(s in 1:2) {
    set.seed(s)
    png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
    png(png_file)
    p1 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label))
    print(p1)
    dev.off()

    random_seq = c(random_seq, rnorm(1))
  }

  # The random numbers are expected to be all different
  expect_true(length(unique(random_seq))==length(random_seq))
})

test_that("calling geom_label_repel without seed creates different plots", {
  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # Make a plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
  png(png_file)
  p1 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label))
  print(p1)
  grid.force()
  pos1 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Make a second plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed2")
  png(png_file)
  p2 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label))
  print(p2)
  grid.force()
  pos2 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Confirm that the label positions are identical.
  expect_true(nrow(pos1) == nrow(dat1))
  expect_true(nrow(pos2) == nrow(dat1))
  expect_true(!identical(pos1, pos2))
})

test_that("calling geom_label_repel with seed creates identical plots", {
  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # Make a plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
  png(png_file)
  p1 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label), seed = 10)
  print(p1)
  grid.force()
  pos1 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Make a second plot with no seed and get the label positions.
  png_file <- withr::local_tempfile(pattern = "testthat_test-seed2")
  png(png_file)
  p2 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label), seed = 10)
  print(p2)
  grid.force()
  pos2 <- pos_df(grid.get("textrepelgrob", grep = TRUE, global = TRUE))
  dev.off()

  # Confirm that the label positions are identical.
  expect_true(nrow(pos1) == nrow(dat1))
  expect_true(nrow(pos2) == nrow(dat1))
  expect_equal(pos1, pos2, tolerance = 1e-3)
})

test_that("calling geom_label_repel without seed does not remove entropy", {
  ix <- seq(1, nrow(mtcars), 4)
  dat1 <- mtcars[ix,]
  dat1$label <- rownames(mtcars)[ix]

  # One random number will be generated after each plot and accumulated
  random_seq = c()
  for(s in 1:2) {
    set.seed(s)
    png_file <- withr::local_tempfile(pattern = "testthat_test-seed1")
    png(png_file)
    p1 <- ggplot(dat1) + geom_label_repel(aes(wt, mpg, label = label))
    print(p1)
    dev.off()

    random_seq = c(random_seq, rnorm(1))
  }

  # The random numbers are expected to be all different
  expect_true(length(unique(random_seq))==length(random_seq))
})

