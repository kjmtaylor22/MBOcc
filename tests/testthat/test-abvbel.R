test_that("vector is binary", {
  expect_in(u.d(seq(0,10), stat="mean"), c(0,1))
  expect_in(u.d(seq(0,10), stat="median"), c(0,1))
  expect_in(u.d(seq(1,10), stat="logmean"), c(0,1))
})


