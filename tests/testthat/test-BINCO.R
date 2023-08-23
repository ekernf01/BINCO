test_that("valley.score like, totally, works", {
  testthat::expect_equal(valley.score(((-50):50)^2)[1], 1/2, tolerance = 0.02)
  testthat::expect_equal(valley.score(((-50):25)^2)[1], 2/3, tolerance = 0.02)
  testthat::expect_equal(valley.score(((-50):10)^2)[1], 5/6, tolerance = 0.02)
  testthat::expect_equal(valley.score(((-50):5)^2)[1], 10/11, tolerance = 0.02)
  testthat::expect_equal(valley.score(((-50):2)^2)[1], 25/26, tolerance = 0.02)
})

test_that("u_shape_test works", {
  testthat::expect_equal(u_shape_test(((-50):50)^2)[1], 0)
  testthat::expect_equal(u_shape_test(((-50):25)^2)[1], 0)
  testthat::expect_equal(u_shape_test(((-50):10)^2)[1], 1)
  testthat::expect_equal(u_shape_test(((-50):5)^2)[1], 1)
  testthat::expect_equal(u_shape_test(((-50):2)^2)[1], 1)
  testthat::expect_equal(u_shape_test(((-50):10)^2, vpr = 1)[1], 0)
  testthat::expect_equal(u_shape_test(((-50):5)^2, vpr = 1)[1], 0)
  testthat::expect_equal(u_shape_test(((-50):2)^2, vpr = 1)[1], 0)
})
