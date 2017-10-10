context("mixopt")

test_that("SQP method without partial QR works for normmix data set",{
  L   <- normmix.data$L
  out <- mixsqp(L,pqrtol = 0,verbose = FALSE)
  expect_equal(normmix.data$w,out$x,tolerance = 1e-4)
})

test_that(paste("SQP method with partial QR closely recovers correct",
                "solution for normmix data set"),{
  L   <- normmix.data$L
  out <- mixsqp(L,pqrtol = 1e-8,verbose = FALSE)
  expect_equal(normmix.data$w,out$x,tolerance = 1e-4)
})

