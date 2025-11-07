test_that("Throws reasonable error when character is included in reads vector", {
  expect_error(predict_detection(c(1000, 10000, 10, "f")))
})

test_that("Throws reasonable error when ci > 1", {
  expect_error(predict_detection(c(1000, 10000, 10), ci = 90))
})

test_that("Throws reasonable error when ci < 0", {
  expect_error(predict_detection(c(1000, 10000, 10), ci = -10))
})

test_that("Throws reasonable error when ci is non-numeric", {
  expect_error(predict_detection(c(1000, 10000, 10), ci = "f"))
})

test_that("Works when single read value is used", {
  expect_equal(class(predict_detection(c(1000))), "list")
  expect_equal(names(predict_detection(c(1000))), c("fit", "lwr", "upr"))
})

test_that("Works when zero value in reads vector", {
  expect_equal(class(predict_detection(c(1000, 10000, 10, 0))), "list")
  expect_equal(
    names(predict_detection(c(1000, 10000, 10, 0))),
    c("fit", "lwr", "upr")
  )
})

test_that("Works when only zero is provided to reads argument", {
  expect_equal(class(predict_detection(c(0))), "list")
  expect_equal(names(predict_detection(c(0))), c("fit", "lwr", "upr"))
})

test_that("returns list of fit, lwr and upr, when default model (DLOOP) is used", {
  expect_equal(class(predict_detection(c(1000, 10000, 10))), "list")
  expect_equal(
    names(predict_detection(c(1000, 10000, 10))),
    c("fit", "lwr", "upr")
  )
})

test_that("returns list of fit, lwr and upr for DLOOP model", {
  expect_equal(class(predict_detection(c(1000, 10000, 10), "DLOOP")), "list")
  expect_equal(
    names(predict_detection(c(1000, 10000, 10), "DLOOP")),
    c("fit", "lwr", "upr")
  )
})

test_that("returns list of fit, lwr and upr for DORY model", {
  expect_equal(class(predict_detection(c(1000, 10000, 10), "DORY")), "list")
  expect_equal(
    names(predict_detection(c(1000, 10000, 10), "DORY")),
    c("fit", "lwr", "upr")
  )
})
