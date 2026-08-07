## Auth database setup — schema creation and local test-admin seeding.
## Standalone: not sourced by app.R yet (no auth logic wired in). Invoked
## explicitly, e.g. from docker-compose.yml before the app starts.
##
## DB_PATH convention: the SQLite file lives wherever DB_PATH points —
## ./local-data/auth.sqlite locally (docker-compose bind mount), overridden
## to the gcsfuse mount path (e.g. /mnt/gcs-auth/auth.sqlite) on Cloud Run.

library(DBI)
library(RSQLite)

source("R/db_utils.R")  # hash_password()

# Create the `users` and `catalog_access_requests` tables if they don't already exist.
create_auth_schema <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)

  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con))

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS users (
      email         TEXT PRIMARY KEY,
      password_hash TEXT NOT NULL,
      role          TEXT NOT NULL DEFAULT 'general' CHECK (role IN ('general','catalog_access','admin')),
      created_at    TEXT DEFAULT (datetime('now'))
    )
  ")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS catalog_access_requests (
      id             INTEGER PRIMARY KEY AUTOINCREMENT,
      user_email     TEXT NOT NULL REFERENCES users(email),
      justification  TEXT,
      status         TEXT NOT NULL DEFAULT 'pending' CHECK (status IN ('pending','approved','denied')),
      requested_at   TEXT DEFAULT (datetime('now')),
      decided_by     TEXT,
      decided_at     TEXT
    )
  ")

  dbExecute(con, "
    CREATE UNIQUE INDEX IF NOT EXISTS one_pending_per_user
      ON catalog_access_requests(user_email) WHERE status = 'pending'
  ")

  invisible(TRUE)
}

# Insert one admin user with a known test password, for local development only.
# Guarded by SEED_TEST_ADMIN=true so this can never run against a real database
# by accident.
seed_test_admin <- function(db_path,
                             email    = "admin@test.local",
                             password = "TestAdmin123!") {
  if (!identical(Sys.getenv("SEED_TEST_ADMIN"), "true")) {
    message("SEED_TEST_ADMIN is not 'true' — skipping test admin seed.")
    return(invisible(FALSE))
  }

  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con))

  dbExecute(
    con,
    "INSERT OR IGNORE INTO users (email, password_hash, role) VALUES (?, ?, 'admin')",
    params = list(email, hash_password(password))
  )

  message(sprintf("Seeded test admin '%s' (password: '%s') — local testing only.", email, password))
  invisible(TRUE)
}
