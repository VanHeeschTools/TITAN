## Catalog access / role-request scaffolding — stubs only, not wired into app.R yet.

# Submit a catalog access request for a user against a given study/role.
# Returns the newly created request record (or its id).
submit_catalog_request <- function(user_id, study_id, role, justification = NULL) {
  stop("submit_catalog_request() not yet implemented")
}

# Fetch all pending (unreviewed) catalog access requests, optionally filtered by study.
# Returns a data.frame of pending requests.
get_pending_requests <- function(study_id = NULL) {
  stop("get_pending_requests() not yet implemented")
}

# Approve a pending request, granting the requested role.
# Returns the updated request record.
approve_request <- function(request_id, reviewer_id, notes = NULL) {
  stop("approve_request() not yet implemented")
}

# Deny a pending request.
# Returns the updated request record.
deny_request <- function(request_id, reviewer_id, notes = NULL) {
  stop("deny_request() not yet implemented")
}

# Check whether a user currently holds a given role (optionally scoped to a study).
# Returns TRUE/FALSE.
user_has_role <- function(user_id, role, study_id = NULL) {
  stop("user_has_role() not yet implemented")
}
