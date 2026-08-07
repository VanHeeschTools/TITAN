## Database utility helpers: password hashing and gcsfuse-tolerant read/write wrappers.
##
## Hashing uses scrypt (already a dependency — see R/db_setup.R and the Dockerfile),
## not sodium/bcrypt: it's the same algorithm shinymanager's check_credentials()
## verifies against, so hashes written here stay compatible with login once that's wired in.

library(DBI)
library(RSQLite)
library(scrypt)

hash_password <- function(password) {
  scrypt::hashPassword(password)
}

verify_password <- function(password, hash) {
  isTRUE(scrypt::verifyPassword(hash, password))
}

# Write with retry — tolerates transient gcsfuse write contention (e.g. a brief
# lock while another writer/reader touches the mounted file). Logs a warning
# and returns FALSE without raising if all retries are exhausted; never crashes
# the app. `.dbExecute` is a testing seam (defaults to DBI::dbExecute).
safe_db_write <- function(db_path, sql, params = NULL, max_retries = 3,
                           .dbExecute = DBI::dbExecute) {
  backoff_s <- c(0.2, 0.5, 1)
  attempt   <- 0

  do_write <- function() {
    con <- dbConnect(RSQLite::SQLite(), db_path)
    on.exit(dbDisconnect(con))
    .dbExecute(con, sql, params = params)
  }

  repeat {
    attempt <- attempt + 1
    result  <- tryCatch(do_write(), error = function(e) e)

    if (!inherits(result, "error")) return(result)

    if (attempt > max_retries) {
      warning(sprintf(
        "safe_db_write: giving up after %d attempt(s) — %s",
        attempt, conditionMessage(result)
      ))
      return(invisible(FALSE))
    }

    Sys.sleep(backoff_s[min(attempt, length(backoff_s))])
  }
}

# Plain read wrapper — no retry: reads aren't subject to the write-lock
# contention safe_db_write guards against.
safe_db_read <- function(db_path, sql, params = NULL) {
  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con))
  dbGetQuery(con, sql, params = params)
}

# Builds the `check_credentials` function shinymanager::secure_server() expects:
# function(user, password) -> list(result, user_info, expired, authorized).
# shinymanager's own check_credentials(db_path) assumes its own `credentials`
# table shape (user/password/is_hashed_password/admin/...); our `users` table
# uses different columns (email, password_hash, role), so we implement the
# contract directly against it instead.
make_check_credentials <- function(db_path) {
  function(user, password) {
    row <- safe_db_read(
      db_path,
      "SELECT email, password_hash, role FROM users WHERE email = ?",
      params = list(user)
    )

    if (nrow(row) == 0) return(list(result = FALSE, user_info = NULL))

    list(
      result     = verify_password(password, row$password_hash[1]),
      user_info  = list(user = row$email[1], role = row$role[1]),
      expired    = FALSE,
      authorized = TRUE
    )
  }
}
