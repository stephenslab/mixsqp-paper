context("mixopt")

test_that("SQP method without partial QR works for normmix data set",{
  L   <- normmix.data$L
  k   <- ncol(L)
  x0  <- rep(1,k)/k
  out <- mixsqp(L,x0,pqrtol = 0,verbose = FALSE)
  expect_equal(normmix.data$w,out$x,tolerance = 1e-4)
})

test_that(paste("SQP method with partial QR closely recovers correct",
                "solution for normmix data set"),{
  L   <- normmix.data$L
  k   <- ncol(L)
  x0  <- rep(1,k)/k
  out <- mixsqp(L,x0,pqrtol = 1e-8,verbose = FALSE)
  expect_equal(normmix.data$w,out$x,tolerance = 1e-4)
})

