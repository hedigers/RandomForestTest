context("test of function hypoRF")
require(stats)

test_that("function hypoRF has correct input parameters", {
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), seed="a"), "seed")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), seed=Inf), "seed")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), seed=-1), "seed")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), probability=T), "probability")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), normalapprox=1), "not logical")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), normalapprox="a"), "not logical")
  expect_error(hypoRF(replicate(3, rnorm(10)),
                      replicate(3, rnorm(10))), "data.frame")
  expect_error(hypoRF(data.frame(replicate(5, rnorm(10))),
                      data.frame(replicate(3, rnorm(10)))), "columns")
  expect_error(hypoRF(data.frame(x=rnorm(10),z=rnorm(10)),
                      data.frame(u=rnorm(10), z=rnorm(10))), "colnames")
  expect_error(hypoRF(data.frame(x=rnorm(10),y=rnorm(10)),
                      data.frame(x=rnorm(10), y=rnorm(10))), "reserved")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), K="a"), "not valid")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), K=Inf), "not valid")
  expect_error(hypoRF(data.frame(replicate(3, rnorm(10))),
                      data.frame(replicate(3, rnorm(10))), K=0), "not valid")
})
