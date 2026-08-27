## Auth database setup — schema creation and local test-admin seeding.
## Depends on hash_password() from R/db_utils.R — callers must source that
## first (not self-sourced here with a relative path: that breaks under
## testthat, which runs tests with a different working directory).
##
## DB_PATH convention: the SQLite file lives wherever DB_PATH points —
## ./local-data/auth.sqlite locally (docker-compose bind mount), overridden
## to the gcsfuse mount path (e.g. /mnt/gcs-auth/auth.sqlite) on Cloud Run.

library(DBI)
library(RSQLite)

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

  # Per-study grants — replaces the old all-or-nothing model where
  # role = 'catalog_access' alone unlocked every study. study_id isn't a
  # DB-level FK (the study catalog lives in data/catalog.yaml, not a table),
  # just the same study_id string used throughout the app (STUDY_CATALOG,
  # catalog_study_list, etc.).
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS study_access (
      user_email  TEXT NOT NULL REFERENCES users(email),
      study_id    TEXT NOT NULL,
      granted_by  TEXT,
      granted_at  TEXT DEFAULT (datetime('now')),
      PRIMARY KEY (user_email, study_id)
    )
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
