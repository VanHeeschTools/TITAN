## Admin: pending catalog access requests — approve/deny.
## Visible only to admins; wired into the navbar in app.R via nav_insert()
## once user_has_role(session, "admin") is known post-login (role isn't
## available before login, so this can't be a static nav_panel).

mod_admin_requests_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header(tags$span(icon("user-shield"), " Pending catalog access requests"), class = "fw-semibold"),
    card_body(DTOutput(ns("table")))
  )
}

mod_admin_requests_server <- function(id, db_path) {
  moduleServer(id, function(input, output, session) {

    refresh_trigger <- reactiveVal(0)
    refresh <- function() refresh_trigger(isolate(refresh_trigger()) + 1)

    pending <- reactive({
      refresh_trigger()
      get_pending_requests(db_path)
    })

    output$table <- renderDT({
      df <- pending()

      if (nrow(df) == 0) {
        return(datatable(
          data.frame(Message = "No pending requests."),
          rownames = FALSE, selection = "none", options = list(dom = "t")
        ))
      }

      row_buttons <- function(rid) {
        as.character(tagList(
          tags$button(
            class = "btn btn-success btn-sm", type = "button",
            onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'})",
                              session$ns("approve_click"), rid),
            "Approve"
          ),
          " ",
          tags$button(
            class = "btn btn-danger btn-sm", type = "button",
            onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'})",
                              session$ns("deny_click"), rid),
            "Deny"
          )
        ))
      }

      display <- data.frame(
        Email         = df$user_email,
        Justification = ifelse(is.na(df$justification), "", df$justification),
        Requested     = df$requested_at,
        Actions       = vapply(df$id, row_buttons, character(1)),
        stringsAsFactors = FALSE
      )

      datatable(
        display,
        escape    = FALSE,
        rownames  = FALSE,
        selection = "none",
        class     = "compact hover",
        options   = list(dom = "t", paging = FALSE, ordering = FALSE)
      )
    })

    observeEvent(input$approve_click, {
      result <- approve_request(db_path, input$approve_click, session$userData$user)
      showNotification(result$message, type = if (result$success) "message" else "error")
      refresh()
    })

    observeEvent(input$deny_click, {
      result <- deny_request(db_path, input$deny_click, session$userData$user)
      showNotification(result$message, type = if (result$success) "message" else "error")
      refresh()
    })
  })
}
