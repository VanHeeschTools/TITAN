-- Placeholder schema for catalog access / role-request feature.
-- Not yet applied anywhere; scaffolding only.

CREATE TABLE IF NOT EXISTS catalog_access_requests (
    id              SERIAL PRIMARY KEY,
    user_id         TEXT        NOT NULL,
    study_id        TEXT        NOT NULL,
    role            TEXT        NOT NULL,
    justification   TEXT,
    status          TEXT        NOT NULL DEFAULT 'pending', -- pending | approved | denied
    requested_at    TIMESTAMPTZ NOT NULL DEFAULT now(),
    reviewed_by     TEXT,
    reviewed_at     TIMESTAMPTZ,
    review_notes    TEXT
);

CREATE INDEX IF NOT EXISTS idx_catalog_access_requests_status
    ON catalog_access_requests (status);

CREATE INDEX IF NOT EXISTS idx_catalog_access_requests_user_id
    ON catalog_access_requests (user_id);
