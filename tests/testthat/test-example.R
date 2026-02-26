test_that("vcocena_example returns valid paths", {
  ex <- vcocena_example()
  expect_true(dir.exists(ex$working_directory))
  expect_true(file.exists(file.path(ex$working_directory, "data", ex$counts_seq)))
  expect_true(file.exists(file.path(ex$working_directory, "data", ex$counts_array)))
  expect_true(file.exists(file.path(ex$working_directory, "sample_info", ex$anno_seq)))
  expect_true(file.exists(file.path(ex$working_directory, "sample_info", ex$anno_array)))
  expect_true(file.exists(file.path(ex$reference_dir, ex$tf)))
})
