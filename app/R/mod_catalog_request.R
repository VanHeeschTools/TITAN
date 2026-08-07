## Catalog access request UI — shown instead of the Catalog tab for users
## without the catalog_access role. This module does NOT check the role
## itself: the call site decides whether to render it at all, e.g.
##   if (user_has_role(session, "catalog_access")) real_catalog_ui() else mod_catalog_request_ui("catalog_request")
## (wiring into the actual Catalog tab is a separate step).

mod_catalog_request_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header(tags$span(icon("lock"), " Catalog access"), class = "fw-semibold"),
    card_body(uiOutput(ns("body")))
  )
}

mod_catalog_request_server <- function(id, db_path) {
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

    output$body <- renderUI({
      status <- request_status()

      if (identical(status, "pending")) {
        tags$p(class = "text-warning mb-0", icon("hourglass-half"),
               " Your catalog access request is pending review.")
      } else if (identical(status, "approved")) {
        tags$p(class = "text-success mb-0", icon("circle-check"),
               " Your catalog access request was approved — log out and back in to pick up your new role.")
      } else {
        tagList(
          if (identical(status, "denied"))
            tags$p(class = "text-danger",
                   "Your previous request was denied. You can submit a new one below."),
          tags$p("Catalog access is required to browse study data. Tell us why you need it:"),
          textAreaInput(session$ns("justification"), label = NULL, rows = 4, width = "100%",
                        placeholder = "e.g. project name, PI, intended use of the catalog"),
          actionButton(session$ns("submit"), "Request Access", class = "btn-primary")
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
  })
}
