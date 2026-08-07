## Catalog access request management: submit / list / approve / deny, and a
## session-level role check. Backed by DBI/SQLite via R/db_utils.R
## (safe_db_write / safe_db_read); schema in R/db_setup.R.

# Submit a pending catalog access request for `email`. Rejects with a friendly
# message (does not throw) if a pending request already exists for that email.
# The pre-check below handles the common case without wasting safe_db_write's
# retries on a guaranteed-to-fail insert; the partial unique index
# (one_pending_per_user ON catalog_access_requests(user_email) WHERE
# status='pending') is the actual correctness guard against a race between
# the check and the insert — if it fires, the write fails and we report the
# same friendly message rather than a generic error.
submit_catalog_request <- function(db_path, email, justification = NULL) {
  existing <- safe_db_read(
    db_path,
    "SELECT 1 FROM catalog_access_requests WHERE user_email = ? AND status = 'pending'",
    params = list(email)
  )
  if (nrow(existing) > 0) {
    return(list(success = FALSE, message = "You already have a pending catalog access request."))
  }

  result <- safe_db_write(
    db_path,
    "INSERT INTO catalog_access_requests (user_email, justification) VALUES (?, ?)",
    params = list(email, justification)
  )

  if (isFALSE(result)) {
    return(list(success = FALSE, message = "You already have a pending catalog access request."))
  }

  list(success = TRUE, message = "Request submitted — an admin will review it shortly.")
}

# All rows with status = 'pending', oldest first.
get_pending_requests <- function(db_path) {
  safe_db_read(
    db_path,
    "SELECT * FROM catalog_access_requests WHERE status = 'pending' ORDER BY requested_at"
  )
}

# Approve a pending request: marks it 'approved' AND grants the requester the
# 'catalog_access' role. Not a cross-table SQLite transaction (not needed
# here) — the two writes run in sequence. If the role grant fails after the
# status update already succeeded, that mismatch is logged loudly (warning +
# the exact fix-up SQL) so it can be corrected manually rather than silently
# leaving the request marked approved without the role actually granted.
approve_request <- function(db_path, request_id, admin_email) {
  pending <- safe_db_read(
    db_path,
    "SELECT user_email FROM catalog_access_requests WHERE id = ? AND status = 'pending'",
    params = list(request_id)
  )
  if (nrow(pending) == 0) {
    return(list(success = FALSE, message = "No pending request with that id."))
  }
  user_email <- pending$user_email[1]

  status_result <- safe_db_write(
    db_path,
    "UPDATE catalog_access_requests SET status = 'approved', decided_by = ?, decided_at = datetime('now') WHERE id = ?",
    params = list(admin_email, request_id)
  )
  if (isFALSE(status_result)) {
    return(list(success = FALSE, message = "Could not update request status — please try again."))
  }

  role_result <- safe_db_write(
    db_path,
    "UPDATE users SET role = 'catalog_access' WHERE email = ?",
    params = list(user_email)
  )
  if (isFALSE(role_result)) {
    warning(sprintf(
      paste0(
        "approve_request: request %s marked 'approved' for '%s' but the role grant FAILED. ",
        "Fix manually: UPDATE users SET role='catalog_access' WHERE email='%s';"
      ),
      request_id, user_email, user_email
    ))
    return(list(
      success = FALSE,
      message = sprintf(
        "Request approved but granting the role failed — contact an admin to fix '%s' manually.",
        user_email
      )
    ))
  }

  list(success = TRUE, message = sprintf("Approved — '%s' now has catalog_access.", user_email))
}

# Deny a pending request. Only updates the request row — the user's role is
# left untouched.
deny_request <- function(db_path, request_id, admin_email) {
  pending <- safe_db_read(
    db_path,
    "SELECT id FROM catalog_access_requests WHERE id = ? AND status = 'pending'",
    params = list(request_id)
  )
  if (nrow(pending) == 0) {
    return(list(success = FALSE, message = "No pending request with that id."))
  }

  result <- safe_db_write(
    db_path,
    "UPDATE catalog_access_requests SET status = 'denied', decided_by = ?, decided_at = datetime('now') WHERE id = ?",
    params = list(admin_email, request_id)
  )
  if (isFALSE(result)) {
    return(list(success = FALSE, message = "Could not deny request — please try again."))
  }

  list(success = TRUE, message = "Request denied.")
}

# Check the role granted at login. Reads session$userData$role, kept in sync
# with the shinymanager auth reactive by an observe() in app.R's server().
user_has_role <- function(session, role) {
  isTRUE(session$userData$role == role)
}
