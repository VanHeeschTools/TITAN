## Formatting utility functions — no Shiny or reactive dependencies.

fmt1     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.1f%%", x)
fmt2     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.3f", x)
bool_fmt <- function(x) if (is.na(x) || length(x) == 0) "—" else if (isTRUE(x)) "Yes ✓" else "No"
