## Peptide matching and protein-sequence display helpers.
## Depends on: stringr (loaded in global.R); html_attr_escape, build_pep_popover defined here.

## Exact substring matching via an integer-encoded k-mer seed index.
## Building a character-keyed index (e.g. via split() on k-mer strings) does not
## scale to tens of millions of k-mers - R's string hashing/sorting is the bottleneck.
## Encoding each k-mer as an integer lets grouping use a fast integer radix sort instead.
match_peptides <- function(peptides, orf_tbl, k = 6L) {
  peptides <- unique(trimws(peptides))
  peptides <- peptides[nchar(peptides) >= 8]
  if (length(peptides) == 0) return(NULL)
  canonical <- c("ORF-annotated", "NC-variant")

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

  # ---- build index: integer k-mer hash -> ORF row indices (one-time cost per call) ----
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
