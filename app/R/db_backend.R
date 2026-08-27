## Backend selection and dialect helpers for the auth DB.
##
## DB_BACKEND=sqlite (default) — local file at db_path, as used everywhere
## today (docker-compose bind mount locally, gcsfuse mount on Cloud Run).
## DB_BACKEND=postgres — networked DB via PG* env vars (PGHOST, PGPORT,
## PGDATABASE, PGUSER, PGPASSWORD); db_path is ignored. On Cloud Run, PGHOST
## is the Cloud SQL Auth Proxy / native connector unix socket path
## (/cloudsql/<connection_name>).
##
## Every call site keeps using safe_db_read()/safe_db_write() with `?`
## placeholders unchanged (translate_placeholders() rewrites them to `$1,
## $2, ...` for Postgres) — this file is the only place that knows which
## backend is active.

library(DBI)

db_backend <- function() {
  backend <- tolower(Sys.getenv("DB_BACKEND", "sqlite"))
  if (!backend %in% c("sqlite", "postgres")) {
    stop(sprintf("Unknown DB_BACKEND '%s' - must be 'sqlite' or 'postgres'.", backend))
  }
  backend
}

get_db_connection <- function(db_path) {
  if (db_backend() == "postgres") {
    library(RPostgres)
    dbConnect(
      RPostgres::Postgres(),
      host     = Sys.getenv("PGHOST"),
      port     = as.integer(Sys.getenv("PGPORT", "5432")),
      dbname   = Sys.getenv("PGDATABASE"),
      user     = Sys.getenv("PGUSER"),
      password = Sys.getenv("PGPASSWORD")
    )
  } else {
    library(RSQLite)
    dbConnect(RSQLite::SQLite(), db_path)
  }
}

# Rewrites sequential `?` placeholders to Postgres `$1, $2, ...`; no-op for
# SQLite (RSQLite accepts `?` natively). None of this app's queries embed a
# literal '?' character in the SQL text itself, so a straightforward
# left-to-right substitution is safe here.
translate_placeholders <- function(sql) {
  if (db_backend() != "postgres") return(sql)
  i <- 0
  while (grepl("?", sql, fixed = TRUE)) {
    i <- i + 1
    sql <- sub("?", paste0("$", i), sql, fixed = TRUE)
  }
  sql
}

# "Insert, or do nothing if it already exists" - dialect differs (SQLite:
# INSERT OR IGNORE; Postgres: ON CONFLICT ... DO NOTHING). All `cols` are
# bound as `?`/`$n` placeholders; pass literal values (e.g. a fixed role) as
# ordinary params rather than embedding them in the SQL text.
sql_insert_or_ignore <- function(table, cols, conflict_cols) {
  col_list     <- paste(cols, collapse = ", ")
  placeholders <- paste(rep("?", length(cols)), collapse = ", ")
  if (db_backend() == "postgres") {
    sprintf(
      "INSERT INTO %s (%s) VALUES (%s) ON CONFLICT (%s) DO NOTHING",
      table, col_list, placeholders, paste(conflict_cols, collapse = ", ")
    )
  } else {
    sprintf("INSERT OR IGNORE INTO %s (%s) VALUES (%s)", table, col_list, placeholders)
  }
}
