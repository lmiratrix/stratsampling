

test_that("get col or vector works", {
  
  dat = tibble( A = 1:5, B = letters[1:5] )
  
  expect_equal( get_col_or_vector( dat, "A" ),
                dat$A )
  
  expect_equal( get_col_or_vector( dat, dat$A ),
                dat$A )
  
  expect_equal( get_col_or_vector( NULL, dat$A ),
                dat$A )
  
  expect_equal( get_col_or_vector( dat, dat$B ),
                dat$B )
  expect_equal( get_col_or_vector( dat, "B" ),
                dat$B )
  expect_equal( get_col_or_vector( NULL, dat$B ),
                dat$B )

  expect_error( get_col_or_vector( dat, "C" ) )  
  #expect_error( get_col_or_vector( NULL, "C" ) )  
  
})
