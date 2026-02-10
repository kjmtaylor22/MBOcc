list1 <- list(X=(1:11), Y=(1:11), Z=(1:11))
list2 <- list(W=(1:10), Y=(1:10), Z=(1:10))
list3 <- list(A=as.data.frame(list1), B=as.data.frame(list1), C=as.data.frame(list1))
list4 <- list(D=as.data.frame(list2), B=as.data.frame(list2), C=as.data.frame(list2))
list5 <- list(K=list3, L=list3, M=list3)
list6 <- list(G=list3, H=list3, I=list4)

# combining lists of vectors
test_that("list of vectors combines", {
  expect_type(recombine(list1, list2, c), "list")
  expect_all_true(unlist(lapply(recombine(list1, list2, c), is.vector)))

  expect_type(recombine(list1, list2, append), "list")
  expect_all_true(unlist(lapply(recombine(list1, list2, append), is.vector)))
})

# combining lists of data frames
test_that("list of data frames combines", {
  expect_type(recombine(list3, list3, cbind), "list")
  expect_all_true(unlist(lapply(recombine(list3, list3, cbind), is.data.frame)))

  expect_type(recombine(list3, list3, rbind), "list")
  expect_all_true(unlist(lapply(recombine(list3, list3, rbind), is.data.frame)))

  expect_error(recombine(list3, list4, rbind))
  expect_error(recombine(list3, list4, cbind))
})

# combining lists of data frames with different lengths and/or columns
test_that("lists of data frames of different lengths combine", {
  expect_all_true(unlist(lapply(recombine(list3, list4, c), is.list)))
  expect_false(all(unlist(lapply(recombine(list3, list4, c), is.data.frame))))

  expect_all_true(unlist(lapply(recombine(list3, list4, append), is.list)))
  expect_false(all(unlist(lapply(recombine(list3, list4, append), is.data.frame))))
})

# using loop method to combine lists in a list
test_that("list of lists is merged into list of data frames",{
  expect_type(recombine(list5, FUN=cbind, recursive=T), "list")

  expect_all_true(unlist(lapply(recombine(list5, FUN=cbind, recursive=T), is.data.frame)))
  expect_all_true(unlist(lapply(recombine(list5, FUN=rbind, recursive=T), is.data.frame)))

  expect_error(recombine(list6, FUN=rbind, recursive=T))
  expect_error(recombine(list6, FUN=cbind, recursive=T))
})
