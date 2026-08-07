source(testthat::test_path("..", "..", "R", "db_utils.R"))

test_that("hash_password / verify_password round-trip", {
  pwd  <- "correct horse battery staple"
  hash <- hash_password(pwd)

  expect_true(verify_password(pwd, hash))
  expect_false(verify_password("wrong password", hash))
})

test_that("hash_password salts — same password hashes differently each call", {
  pwd <- "same password"
  expect_false(identical(hash_password(pwd), hash_password(pwd)))
})

test_that("safe_db_write retries on transient failure and returns the eventual result", {
  n_calls <- 0
  flaky_execute <- function(...) {
    n_calls <<- n_calls + 1
    if (n_calls <= 2) stop("database is locked")
    42L
  }

  db_path <- tempfile(fileext = ".sqlite")
  DBI::dbDisconnect(DBI::dbConnect(RSQLite::SQLite(), db_path))
  on.exit(unlink(db_path))

  result <- safe_db_write(
    db_path, "INSERT INTO x VALUES (1)",
    max_retries = 3, .dbExecute = flaky_execute
  )

  expect_equal(result, 42L)
  expect_equal(n_calls, 3)
})

test_that("safe_db_write gives up, warns, and does not error after exhausting retries", {
  always_fails <- function(...) stop("database is locked")

  db_path <- tempfile(fileext = ".sqlite")
  DBI::dbDisconnect(DBI::dbConnect(RSQLite::SQLite(), db_path))
  on.exit(unlink(db_path))

  result <- NULL
  expect_warning(
    result <- safe_db_write(
      db_path, "INSERT INTO x VALUES (1)",
      max_retries = 2, .dbExecute = always_fails
    ),
    "safe_db_write"
  )

  expect_false(result)
})
