## Peptide matching and protein-sequence display helpers.
## Depends on: stringr (loaded in global.R); html_attr_escape, build_pep_popover defined here.

## Build the ORF-side half of the k-mer seed index: integer k-mer hash -> sorted
## group of ORF row indices whose protein_seq contains that k-mer. This is the
## part of match_peptides() that only depends on orf_tbl (not on the peptides
## being searched), so for a large static reference table - e.g. the Gencode
## cross-match table, which never changes between peptide uploads - it can be
## built once (scripts/reference_prep/03_prep_gencode_kmer_index.R) and reused
## via match_peptides(..., index = )  instead of being rebuilt on every call.
build_pep_kmer_index <- function(orf_tbl, k = 6L) {
  seqs <- orf_tbl$protein_seq
  base <- 27L   # 26 letters + '*' (stop codon char, if present in protein_seq)

  code_lookup <- integer(256)
  code_lookup[utf8ToInt(paste0(LETTERS, collapse = ""))] <- 0:25
  code_lookup[utf8ToInt("*")] <- 26L

  # integer hash of every k-mer in a sequence via a Horner rolling scheme
  seq_kmer_hashes <- function(s, k) {
    codes <- code_lookup[utf8ToInt(s)]
    n <- length(codes)
    m <- n - k + 1L
    if (m < 1L) return(integer(0))
    h <- integer(m)
    for (j in 0:(k - 1L)) h <- h * base + codes[(1L + j):(m + j)]
    h
  }

  per_orf  <- lapply(seq_along(seqs), function(i) seq_kmer_hashes(seqs[i], k))
  n_kmers  <- lengths(per_orf)
  all_hash <- unlist(per_orf, use.names = FALSE)
  orf_rep  <- rep(seq_along(seqs), n_kmers)

  ord         <- order(all_hash)          # integer radix sort - fast
  hash_sorted <- all_hash[ord]
  orf_sorted  <- orf_rep[ord]
  grp_start_l <- c(TRUE, hash_sorted[-1] != hash_sorted[-length(hash_sorted)])
  grp_id      <- cumsum(grp_start_l)
  uniq_hash   <- hash_sorted[grp_start_l]
  grp_end     <- cumsum(tabulate(grp_id))
  grp_begin   <- c(1L, head(grp_end, -1) + 1L)

  list(k = k, base = base, orf_ids = orf_tbl$orf_id,
       uniq_hash = uniq_hash, grp_begin = grp_begin, grp_end = grp_end,
       orf_sorted = orf_sorted)
}

## Exact substring matching via an integer-encoded k-mer seed index.
## Building a character-keyed index (e.g. via split() on k-mer strings) does not
## scale to tens of millions of k-mers - R's string hashing/sorting is the bottleneck.
## Encoding each k-mer as an integer lets grouping use a fast integer radix sort instead.
##
## `index`, if supplied, must come from build_pep_kmer_index(orf_tbl, k). Reusing a
## prebuilt index skips rebuilding it from orf_tbl$protein_seq on every call - the
## dominant cost when orf_tbl is large and static (e.g. the Gencode cross-match
## table). Falls back to building fresh (with a warning) if index$orf_ids doesn't
## match orf_tbl$orf_id (stale index).
match_peptides <- function(peptides, orf_tbl, k = 6L, index = NULL) {
  peptides <- unique(trimws(peptides))
  peptides <- peptides[nchar(peptides) >= 8]
  if (length(peptides) == 0) return(NULL)
  canonical <- c("ORF-annotated", "NC-variant")

  if (!is.null(index) && !identical(index$orf_ids, orf_tbl$orf_id)) {
    warning("match_peptides(): prebuilt index doesn't match orf_tbl (stale? rerun 03_prep_gencode_kmer_index.R) — rebuilding.")
    index <- NULL
  }
  if (is.null(index)) index <- build_pep_kmer_index(orf_tbl, k = k)

  seqs      <- orf_tbl$protein_seq
  base      <- index$base
  k         <- index$k
  uniq_hash <- index$uniq_hash
  grp_begin <- index$grp_begin
  grp_end   <- index$grp_end
  orf_sorted <- index$orf_sorted

  # ---- vectorised peptide-side seed hashing (no per-peptide R calls) ----
  char_to_code <- function(chars) {
    code <- match(chars, LETTERS) - 1L
    code[chars == "*"] <- 26L
    code
  }
  seed_hash <- integer(length(peptides))
  for (j in seq_len(k)) seed_hash <- seed_hash * base + char_to_code(substr(peptides, j, j))

  gi    <- match(seed_hash, uniq_hash)   # NA where seed never occurs anywhere
  valid <- which(!is.na(gi))
  if (length(valid) == 0) return(NULL)

  g_valid <- gi[valid]
  len     <- grp_end[g_valid] - grp_begin[g_valid] + 1L

  # ---- build ALL (peptide, candidate ORF) pairs at once, instead of looping per peptide ----
  pep_idx_rep <- rep(valid, len)
  offsets     <- unlist(lapply(seq_along(valid), function(i)
                   grp_begin[g_valid[i]]:grp_end[g_valid[i]]), use.names = FALSE)
  cand_orf    <- orf_sorted[offsets]

  # ---- ONE vectorised verification pass instead of one grepl() call per peptide ----
  keep        <- stringi::stri_detect_fixed(seqs[cand_orf], peptides[pep_idx_rep])
  pep_idx_rep <- pep_idx_rep[keep]
  cand_orf    <- cand_orf[keep]
  if (length(cand_orf) == 0) return(NULL)

  # de-duplicate (peptide, orf) pairs that can arise from overlapping seed hits
  combo_key   <- as.double(pep_idx_rep) * (length(seqs) + 1) + cand_orf
  dedup       <- !duplicated(combo_key)
  pep_idx_rep <- pep_idx_rep[dedup]
  cand_orf    <- cand_orf[dedup]

  # ---- ONE bulk subset instead of one per peptide ----
  matched <- orf_tbl[cand_orf, , drop = FALSE]
  matched$matched_peptide <- peptides[pep_idx_rep]

  # ---- vectorised canonical-biotype filter (peptides matching a canonical
  # biotype are not evidence for ncORFs) ----
  is_canon      <- matched$orf_biotype_single %in% canonical
  pep_has_canon <- ave(is_canon, matched$matched_peptide, FUN = any)
  matched[!pep_has_canon | is_canon, , drop = FALSE]
}

pill_badge <- function(text, color = "primary") {
  tags$span(class = paste0("badge rounded-pill bg-", color, " me-1"), text)
}

orf_id_labels <- function(tbl) {
  paste0(tbl$gene_name, "_", tbl$orf_biotype_single, "_",
         tbl$protein_length, "aa_",
         tbl$chr, ":", tbl$orf_start, "-", tbl$orf_end, "_",
         tbl$start_codon)
}

# Escape a string for use inside an HTML attribute value (data-bs-content etc.)
html_attr_escape <- function(s) {
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub('"',  "&quot;", s, fixed = TRUE)
  s <- gsub("'",  "&#39;",  s, fixed = TRUE)
  s
}

# Build an HTML table for the Bootstrap popover from a data.frame of MS rows
build_pep_popover <- function(rows_df) {
  if (is.null(rows_df) || nrow(rows_df) == 0L || ncol(rows_df) == 0L) return("")
  make_table <- function(row) {
    cells <- paste(vapply(colnames(row), function(col) {
      sprintf('<tr><td class="pep-tt-key">%s</td><td class="pep-tt-val">%s</td></tr>',
              col, row[[col]])
    }, character(1)), collapse = "")
    sprintf('<table class="pep-tt-table">%s</table>', cells)
  }
  parts <- vapply(seq_len(nrow(rows_df)), function(i) make_table(rows_df[i, , drop = FALSE]),
                  character(1))
  paste(parts, collapse = "<hr class='my-1'>")
}

# Render protein sequence as HTML with per-peptide colour highlights,
# alignment rows, and Bootstrap popover tooltips showing MS data on hover.
# pep_info: named list  peptide → data.frame of MS rows (for popover)
render_protein_seq_html <- function(seq, pep_list, pep_info = list()) {
  PEP_COLS <- c("#2F3D46", "#D4850A", "#8E44AD", "#C0392B", "#0097A7")
  pep_list <- unique(pep_list[!is.na(pep_list) & nzchar(pep_list)])
  n_chars  <- nchar(seq)
  seq_v    <- strsplit(seq, "")[[1]]

  # Mark which peptide (1-indexed) first covers each position
  coverage   <- integer(n_chars)
  pep_starts <- vector("list", length(pep_list))
  for (pi in seq_along(pep_list)) {
    m <- gregexpr(pep_list[[pi]], seq, fixed = TRUE)[[1]]
    if (m[1L] > 0L) {
      pep_starts[[pi]] <- m
      plen <- nchar(pep_list[[pi]])
      for (s in m) for (j in s:min(s + plen - 1L, n_chars)) if (!coverage[j]) coverage[j] <- pi
    }
  }

  BLOCK  <- 60L
  INDENT <- "     "   # 4-digit line number + 1 space

  blocks <- vapply(seq_len(ceiling(n_chars / BLOCK)), function(b) {
    i0 <- (b - 1L) * BLOCK + 1L
    i1 <- min(b * BLOCK, n_chars)

    # Sequence row - use RLE to group runs of same peptide into one span
    cov_range <- coverage[i0:i1]
    seq_range <- seq_v[i0:i1]
    r_len  <- rle(cov_range)$lengths
    r_val  <- rle(cov_range)$values
    r_end  <- cumsum(r_len)
    r_start <- c(1L, r_end[-length(r_end)] + 1L)

    seq_parts <- mapply(function(rs, re, cv) {
      chars <- paste(seq_range[rs:re], collapse = "")
      if (cv == 0L) return(chars)
      col <- PEP_COLS[(cv - 1L) %% length(PEP_COLS) + 1L]
      pep <- pep_list[[cv]]
      info      <- pep_info[[pep]]
      n_records <- if (!is.null(info) && nrow(info) > 0L) nrow(info) else 0L
      tt_title   <- html_attr_escape(sprintf("MS data (%d record%s)", n_records, if (n_records == 1L) "" else "s"))
      tt_content <- html_attr_escape(build_pep_popover(info))
      sprintf(
        '<span class="pep-hit" style="background:%s33;color:%s;font-weight:bold;" data-bs-toggle="popover" data-bs-html="true" data-bs-placement="top" data-bs-trigger="hover focus" data-bs-title="%s" data-bs-content="%s">%s</span>',
        col, col, tt_title, tt_content, chars
      )
    }, r_start, r_end, r_val, SIMPLIFY = TRUE)
    seq_line <- sprintf('<span class="seq-pos">%4d</span> %s', i0, paste(seq_parts, collapse = ""))

    # One alignment row per peptide overlapping this block
    aln_rows <- vapply(seq_along(pep_list), function(pi) {
      starts <- pep_starts[[pi]]
      if (is.null(starts) || starts[1L] < 0L) return("")
      pv   <- strsplit(pep_list[[pi]], "")[[1]]
      col  <- PEP_COLS[(pi - 1L) %% length(PEP_COLS) + 1L]
      aln  <- rep(" ", i1 - i0 + 1L)
      for (s in starts) {
        if (s + length(pv) - 1L < i0 || s > i1) next
        for (j in seq_along(pv)) {
          ap <- s + j - 1L
          if (ap >= i0 && ap <= i1) aln[ap - i0 + 1L] <- pv[j]
        }
      }
      if (all(aln == " ")) return("")
      spans <- vapply(aln, function(ch) {
        if (ch == " ") return(ch)
        sprintf('<span style="color:%s;font-weight:bold;">%s</span>', col, ch)
      }, character(1))
      paste0(INDENT, paste(spans, collapse = ""))
    }, character(1))
    aln_rows <- aln_rows[nzchar(aln_rows)]

    paste(c(seq_line, aln_rows, ""), collapse = "\n")
  }, character(1))

  legend_html <- if (length(pep_list)) {
    badges <- paste(vapply(seq_along(pep_list), function(pi) {
      col <- PEP_COLS[(pi - 1L) %% length(PEP_COLS) + 1L]
      sprintf('<code class="seq-legend-badge" style="background:%s22;color:%s;border:1px solid %s55;">%s</code>',
              col, col, col, pep_list[[pi]])
    }, character(1)), collapse = " ")
    sprintf('<div class="seq-legend"><span class="fw-semibold text-muted small me-2">%d MS peptide%s; hover to see MS data:</span>%s</div>',
            length(pep_list), if (length(pep_list) > 1L) "s" else "", badges)
  } else {
    '<p class="text-muted small mb-1">No MS peptides identified for this ORF.</p>'
  }

  HTML(paste0(legend_html,
              '<div class="titan-protein-seq">',
              paste(blocks, collapse = "\n"),
              '</div>'))
}

## Depends on: match_peptides() (this file); gencode_orf_tbl, gencode_kmer_index
## (global.R, static reference objects, safe to reference directly).
##
## Runs the full peptide<->ORF matching pipeline shared by the manual EXPLORE
## TARGETS click and the automatic bundled-peptide matching path (Study
## Library studies without precomputed_pep_hits): in-house match_peptides()
## (unless `hits` is pre-supplied - the precomputed-hits path), ORF-annotated-
## priority collapse, Gencode cross-match, then join to `meta` (one row per
## matched_peptide, MS metadata). Previously this ~80-line body was
## copy-pasted verbatim between the two call sites in app.R; extracted here so
## a third call site (auto-matching when no precomputed hits exist) didn't
## mean a third copy.
##
## `progress`, if given, is called as progress(value, detail) at each stage -
## callers wrap this in their own withProgress() and pass e.g.
## function(v, d) setProgress(v, detail = d).
## `canon_hits`, if supplied, is dat$precomputed_canonical_pep_hits (from
## prepare_titan_inputs.R's paths.reference_proteome) - the caller decides
## whether it's valid for this peptide set (same fingerprint check already
## used for `hits`) and passes NULL when it isn't.
build_pep_orf_matches <- function(peptides, meta, orf_tbl, hits = NULL,
                                   progress = function(value, detail) NULL,
                                   canon_hits = NULL) {
  if (is.null(hits)) {
    progress(0.3, "Matching peptides to ORFs…")
    hits <- match_peptides(peptides, orf_tbl)
    # For peptides that match an ORF-annotated (canonical) entry AND non-canonical
    # entries, discard the non-canonical hits — they are canonical peptide evidence.
    if (!is.null(hits) && nrow(hits) > 0L) {
      hits <- hits %>%
        group_by(matched_peptide) %>%
        filter(
          if (any(orf_biotype_single == "ORF-annotated", na.rm = TRUE))
            orf_biotype_single == "ORF-annotated"
          else
            TRUE
        ) %>%
        ungroup()
    }
  }

  # Gencode cross-match: annotate in-house hits and add Gencode-only rows
  progress(0.6, "Running Gencode cross-match…")
  if (!is.null(gencode_orf_tbl)) {
    gc_all <- match_peptides(peptides, gencode_orf_tbl, index = gencode_kmer_index)
    if (!is.null(gc_all) && nrow(gc_all) > 0L) {
      # Case (a): build per-peptide summary of matching Gencode ORFs
      gc_summary <- gc_all %>%
        group_by(matched_peptide) %>%
        summarise(
          gencode_match_ids = paste(
            sprintf("%s (%s, %s)", orf_id, gene_name, orf_biotype_single),
            collapse = "; "
          ),
          .groups = "drop"
        )
      if (!is.null(hits) && nrow(hits) > 0L) {
        hits <- hits %>%
          left_join(gc_summary, by = "matched_peptide") %>%
          mutate(gencode_match_ids = replace_na(gencode_match_ids, ""),
                 gencode_only       = FALSE)
        # Case (b): peptides with Gencode hits but no in-house hit → new rows
        gc_only_peps <- setdiff(unique(gc_all$matched_peptide), unique(hits$matched_peptide))
      } else {
        gc_only_peps <- unique(gc_all$matched_peptide)
      }
      if (length(gc_only_peps) > 0L) {
        gc_only_rows <- gc_all %>%
          filter(matched_peptide %in% gc_only_peps) %>%
          mutate(gencode_match_ids = "", gencode_only = TRUE)
        # Populate gene-level expression / GTEx / TCGA metrics by borrowing from
        # any in-house ORF of the same gene (these columns are gene-level, not ORF-level)
        expr_cols <- intersect(
          c("target_expression_num_samples", "target_expression_pct_samples",
            "target_expression_median_TPM", "target_expression_max_TPM",
            "GTEX_max_median_TPM", "GTEX_median_TPM", "GTEX_DE_sig_in_all",
            "GTEX_tumor_only", "GTEX_tumor_enriched", "GTEX_tissues_q3_gt1",
            "TCGA_tumor_num_samples", "TCGA_tumor_pct_samples",
            "TCGA_tumor_median_TPM", "TCGA_tumor_max_TPM",
            "TCGA_normal_num_samples", "TCGA_normal_pct_samples",
            "TCGA_normal_median_TPM", "TCGA_normal_max_TPM"),
          colnames(orf_tbl)
        )
        gene_expr <- orf_tbl %>%
          group_by(gene_id_clean) %>%
          summarise(across(all_of(expr_cols), first), .groups = "drop")
        gc_only_rows <- gc_only_rows %>%
          left_join(gene_expr, by = "gene_id_clean")
        hits <- bind_rows(hits, gc_only_rows)
      }
    } else if (!is.null(hits) && nrow(hits) > 0L) {
      hits <- hits %>% mutate(gencode_match_ids = "", gencode_only = FALSE)
    }
  } else if (!is.null(hits) && nrow(hits) > 0L) {
    hits <- hits %>% mutate(gencode_match_ids = "", gencode_only = FALSE)
  }

  # Reference-proteome cross-match (canonical-only genes with no candidate
  # ORF at all): only for peptides that still have no match at this point.
  if (!is.null(canon_hits)) {
    residual_peps <- if (!is.null(hits) && nrow(hits) > 0L)
      setdiff(peptides, unique(hits$matched_peptide)) else peptides
    canon_rows <- build_canonical_only_rows(residual_peps, canon_hits)
    if (!is.null(canon_rows) && nrow(canon_rows) > 0L)
      hits <- bind_rows(hits, canon_rows)
  }

  progress(0.9, "Joining MS metadata…")
  if (!is.null(hits) && nrow(hits) > 0L) {
    left_join(hits, meta, by = "matched_peptide")
  } else {
    data.frame(orf_id = character(0), matched_peptide = character(0))
  }
}

## Builds ORF-shaped rows for peptides that only matched a canonical reference
## protein (dat$precomputed_canonical_pep_hits, from prepare_titan_inputs.R's
## paths.reference_proteome) — the gene has no candidate ORF at all, so there
## is no ribo-seq/translation evidence to show; target_translation_*/
## ribocrypt_* columns are simply absent here and come back as NA once
## bind_rows() merges these into the main hits table.
## `residual_peps` should already be peptides not covered by any in-house or
## Gencode-cross-match hit — canon_hits itself may still contain other peptides.
build_canonical_only_rows <- function(residual_peps, canon_hits) {
  if (length(residual_peps) == 0L || is.null(canon_hits) || nrow(canon_hits) == 0L) return(NULL)
  rows <- canon_hits[canon_hits$matched_peptide %in% residual_peps, , drop = FALSE]
  if (nrow(rows) == 0L) return(NULL)
  rows$orf_id             <- paste0("CANON_", rows$gene_id)
  rows$gene_name          <- rows$gene_symbol
  rows$orf_biotype_single <- "Canonical (no ORF candidate)"
  rows$gencode_match_ids  <- ""
  rows$gencode_only       <- FALSE
  rows$canonical_only     <- TRUE
  # No candidate ORF exists for these genes, so there's no ORF-level coordinate
  # metadata at all — protein_seq (and protein_length derived from it) IS
  # available though, straight from the reference proteome match, so the ORF
  # Detail view can still show the canonical sequence. Set explicitly rather
  # than relying on bind_rows() to backfill NA — if `hits` is itself NULL/empty
  # (a peptide set with zero in-house or Gencode matches), there's no other
  # side of the bind for these columns to inherit, and prioritised_data()
  # requires them.
  rows$protein_length <- nchar(rows$protein_seq)
  for (col in c("chr", "orf_start", "orf_end", "strand", "start_codon", "gene_biotype"))
    rows[[col]] <- NA
  distinct(rows, orf_id, matched_peptide, .keep_all = TRUE)
}
