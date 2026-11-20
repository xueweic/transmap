test_that("package metadata is available", {
  expect_s3_class(utils::packageVersion("transmap"), "package_version")
})
