source(testthat::test_path("..", "..", "R", "db_utils.R"))
source(testthat::test_path("..", "..", "R", "db_setup.R"))
source(testthat::test_path("..", "..", "R", "catalog_access.R"))

make_test_db <- function(email = "alice@example.com") {
  db_path <- tempfile(fileext = ".sqlite")
  create_auth_schema(db_path)
  safe_db_write(
    db_path,
    "INSERT INTO users (email, password_hash, role) VALUES (?, ?, 'general')",
    params = list(email, "dummyhash")
  )
  db_path
}

test_that("submit_catalog_request rejects a duplicate pending request", {
  db_path <- make_test_db()
  on.exit(unlink(db_path))

  first  <- submit_catalog_request(db_path, "alice@example.com", "need it for X")
  second <- submit_catalog_request(db_path, "alice@example.com", "need it again")

  expect_true(first$success)
  expect_false(second$success)

  rows <- safe_db_read(
    db_path, "SELECT * FROM catalog_access_requests WHERE user_email = ?",
    params = list("alice@example.com")
  )
  expect_equal(nrow(rows), 1)
})

test_that("approve_request updates both the request status and the user's role", {
  db_path <- make_test_db()
  on.exit(unlink(db_path))

  submit_catalog_request(db_path, "alice@example.com", "need it")
  req_id <- get_pending_requests(db_path)$id[1]

  result <- approve_request(db_path, req_id, "admin@example.com")
  expect_true(result$success)

  req <- safe_db_read(
    db_path, "SELECT status, decided_by FROM catalog_access_requests WHERE id = ?",
    params = list(req_id)
  )
  expect_equal(req$status[1], "approved")
  expect_equal(req$decided_by[1], "admin@example.com")

  user <- safe_db_read(db_path, "SELECT role FROM users WHERE email = ?", params = list("alice@example.com"))
  expect_equal(user$role[1], "catalog_access")
})

test_that("deny_request updates only the request, leaving the user's role unchanged", {
  db_path <- make_test_db()
  on.exit(unlink(db_path))

  submit_catalog_request(db_path, "alice@example.com", "need it")
  req_id <- get_pending_requests(db_path)$id[1]

  result <- deny_request(db_path, req_id, "admin@example.com")
  expect_true(result$success)

  req <- safe_db_read(db_path, "SELECT status FROM catalog_access_requests WHERE id = ?", params = list(req_id))
  expect_equal(req$status[1], "denied")

  user <- safe_db_read(db_path, "SELECT role FROM users WHERE email = ?", params = list("alice@example.com"))
  expect_equal(user$role[1], "general")
})
