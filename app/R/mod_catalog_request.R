## Catalog access request UI — shown in place of the study library inside the
## ORF candidates card (see catalog_tab_ui() / orf_source_ui in app.R) for
## users without the catalog_access role. This module does NOT check the
## role itself; the call site decides whether to render it.
##
## No outer card() here: it's embedded inside an existing card_body(), so
## wrapping it in another card would nest cards inside a card.

mod_catalog_request_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("body"))
}

# `on_switch_to_upload`: optional zero-arg callback fired when the user
# clicks "upload your own data instead" — wired in app.R to show_upload_rv(TRUE).
mod_catalog_request_server <- function(id, db_path, on_switch_to_upload = NULL) {
  moduleServer(id, function(input, output, session) {

    user_email <- reactive(session$userData$user)

    # Most recent request status for this user: NULL / "pending" / "approved" / "denied".
    request_status <- reactiveVal(NULL)

    refresh_status <- function() {
      email <- user_email()
      if (is.null(email) || !nzchar(email)) {
        request_status(NULL)
        return(invisible(NULL))
      }
      row <- safe_db_read(
        db_path,
        "SELECT status FROM catalog_access_requests WHERE user_email = ? ORDER BY requested_at DESC LIMIT 1",
        params = list(email)
      )
      request_status(if (nrow(row) == 0) NULL else row$status[1])
    }

    # Load current status once the user's email is available (i.e. on login).
    observeEvent(user_email(), refresh_status(), ignoreNULL = FALSE)

    upload_prompt <- tags$p(
      class = "text-muted small mt-2 mb-0",
      "Don't want to wait? ",
      actionLink(session$ns("switch_to_upload"), "Upload your own data instead.")
    )

    output$body <- renderUI({
      status <- request_status()

      if (identical(status, "pending")) {
        tagList(
          tags$p(class = "text-warning mb-1", icon("hourglass-half"),
                 " Your catalog access request is pending review."),
          upload_prompt
        )
      } else if (identical(status, "approved")) {
        tags$p(class = "text-success mb-0", icon("circle-check"),
               " Your catalog access request was approved — log out and back in to pick up your new role.")
      } else {
        tagList(
          if (identical(status, "denied"))
            tags$p(class = "text-danger",
                   "Your previous request was denied. You can submit a new one below."),
          tags$p(class = "small", "Catalog access is required to browse the study library. Tell us why you need it:"),
          textAreaInput(session$ns("justification"), label = NULL, rows = 3, width = "100%",
                        placeholder = "e.g. project name, PI, intended use of the catalog"),
          actionButton(session$ns("submit"), "Request Access", class = "btn-primary btn-sm"),
          upload_prompt
        )
      }
    })

    observeEvent(input$submit, {
      email <- user_email()
      if (is.null(email) || !nzchar(email)) {
        showNotification("You must be logged in to request access.", type = "error")
        return(invisible(NULL))
      }

      result <- submit_catalog_request(db_path, email, input$justification)

      if (result$success) {
        showNotification(result$message, type = "message", duration = 6)
        refresh_status()
      } else {
        showNotification(result$message, type = "error")
      }
    })

    observeEvent(input$switch_to_upload, {
      if (is.function(on_switch_to_upload)) on_switch_to_upload()
    })
  })
}
