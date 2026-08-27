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

# `email`'s own most recent pending request, if any (0-row data.frame if
# none) - used to show "request pending since ..." instead of re-offering
# the request form (submit_catalog_request() already blocks a second pending
# request server-side; this is what lets the UI show that state up front).
get_own_pending_request <- function(db_path, email) {
  safe_db_read(
    db_path,
    "SELECT * FROM catalog_access_requests WHERE user_email = ? AND status = 'pending' ORDER BY requested_at DESC LIMIT 1",
    params = list(email)
  )
}

# Grant `user_email` access to each study in `study_ids` (INSERT OR IGNORE,
# so re-granting an already-held study is a no-op rather than an error).
grant_study_access <- function(db_path, user_email, study_ids, granted_by) {
  study_ids <- unique(study_ids)
  if (length(study_ids) == 0) {
    return(list(success = TRUE, message = "No studies selected."))
  }
  all_ok <- TRUE
  for (sid in study_ids) {
    result <- safe_db_write(
      db_path,
      "INSERT OR IGNORE INTO study_access (user_email, study_id, granted_by) VALUES (?, ?, ?)",
      params = list(user_email, sid, granted_by)
    )
    if (isFALSE(result)) all_ok <- FALSE
  }
  if (!all_ok) {
    warning(sprintf("grant_study_access: one or more grants failed for '%s' (studies: %s).",
                     user_email, paste(study_ids, collapse = " ")))
    return(list(success = FALSE, message = "Some studies could not be granted — please retry."))
  }
  list(success = TRUE, message = sprintf(
    "Granted access to %d stud%s.", length(study_ids), if (length(study_ids) == 1) "y" else "ies"
  ))
}

# study_id values `user_email` currently has access to (character vector,
# possibly empty).
get_user_study_access <- function(db_path, user_email) {
  safe_db_read(db_path, "SELECT study_id FROM study_access WHERE user_email = ?",
               params = list(user_email))$study_id
}

# Per-study access check for the Study Library's Load buttons - admins bypass
# (same convention as user_has_role()'s admin-bypass), everyone else needs an
# explicit study_access row for this exact study_id.
user_has_study_access <- function(session, db_path, study_id) {
  if (isTRUE(session$userData$role == "admin")) return(TRUE)
  user_email <- session$userData$user
  if (is.null(user_email)) return(FALSE)
  study_id %in% get_user_study_access(db_path, user_email)
}

# Approve a pending request: marks it 'approved', grants the requester access
# to exactly the studies the admin selected (study_ids - picked at approval
# time, not "all studies"), and sets role = 'catalog_access' as a coarse "has
# been approved for something" flag (gates whether the Study Library tab/menu
# item appears at all elsewhere in the app; per-study Load buttons check
# study_access directly via user_has_study_access(), not this role). Not a
# cross-table SQLite transaction (not needed here) — the writes run in
# sequence; any failure after the status update already succeeded is logged
# loudly (warning + the exact fix-up SQL) so it can be corrected manually
# rather than silently leaving a mismatched state.
approve_request <- function(db_path, request_id, admin_email, study_ids) {
  pending <- safe_db_read(
    db_path,
    "SELECT user_email FROM catalog_access_requests WHERE id = ? AND status = 'pending'",
    params = list(request_id)
  )
  if (nrow(pending) == 0) {
    return(list(success = FALSE, message = "No pending request with that id."))
  }
  if (length(study_ids) == 0) {
    return(list(success = FALSE, message = "Select at least one study to grant."))
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

  grant_result <- grant_study_access(db_path, user_email, study_ids, admin_email)
  if (!grant_result$success) {
    warning(sprintf(
      "approve_request: request %s marked 'approved' for '%s' but granting studies FAILED (%s).",
      request_id, user_email, paste(study_ids, collapse = ", ")
    ))
    return(list(success = FALSE,
                message = paste("Request approved but granting access failed:", grant_result$message)))
  }

  # Never downgrades an existing admin - defensive, shouldn't normally be
  # reachable since admins don't file catalog access requests.
  role_result <- safe_db_write(
    db_path,
    "UPDATE users SET role = 'catalog_access' WHERE email = ? AND role != 'admin'",
    params = list(user_email)
  )
  if (isFALSE(role_result)) {
    warning(sprintf(
      paste0(
        "approve_request: studies granted for '%s' but the role update FAILED. ",
        "Fix manually: UPDATE users SET role='catalog_access' WHERE email='%s';"
      ),
      user_email, user_email
    ))
  }

  list(success = TRUE, message = sprintf(
    "Approved — '%s' now has access to %d stud%s.",
    user_email, length(study_ids), if (length(study_ids) == 1) "y" else "ies"
  ))
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
# `role` is a single, mutually-exclusive column (general/catalog_access/admin
# — see the CHECK constraint in db_setup.R), so an admin account can't also
# literally hold 'catalog_access'. Admins bypass every role gate instead,
# matching the normal expectation that admin implies full access.
user_has_role <- function(session, role) {
  actual <- session$userData$role
  isTRUE(actual == role) || isTRUE(actual == "admin")
}
