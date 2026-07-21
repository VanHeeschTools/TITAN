## Overview plotly chart functions.
## Depends on: BIOTYPE_COLORS (global.R), plotly (loaded in global.R).

biotype_bar <- function(df) {
  if (nrow(df) == 0) {
    return(
      plot_ly() %>%
        layout(
          annotations = list(list(text = "No data", x = 0.5, y = 0.5,
                                  xref = "paper", yref = "paper",
                                  showarrow = FALSE, font = list(size = 13, color = "#888"))),
          xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
          paper_bgcolor = "white", plot_bgcolor = "white"
        ) %>% config(displayModeBar = FALSE)
    )
  }
  d <- df %>%
    count(orf_biotype_single, name = "n") %>%
    arrange(desc(n)) %>%
    mutate(orf_biotype_single = factor(orf_biotype_single, levels = rev(orf_biotype_single)))
  cols <- coalesce(BIOTYPE_COLORS[as.character(d$orf_biotype_single)], "#95A5A6")
  plot_ly(d, x = ~n, y = ~orf_biotype_single, type = "bar", orientation = "h",
          marker = list(color = cols),
          text  = ~paste0(n, " (", round(100 * n / sum(n), 1), "%)"),
          textposition = "outside",
          cliponaxis = FALSE,
          hovertemplate = "%{y}: %{x} ORFs<extra></extra>") %>%
    layout(
      xaxis  = list(title = "Number of ORFs", showgrid = TRUE, gridcolor = "#EEF2F7",
                    autorange = TRUE),
      yaxis  = list(title = "", automargin = TRUE),
      paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(l = 5, r = 110, t = 10, b = 30),
      font   = list(family = "Inter", size = 12)
    ) %>% config(displayModeBar = FALSE)
}
