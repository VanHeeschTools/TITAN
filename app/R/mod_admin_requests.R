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

    # Set by the Approve click, read by the modal's own "Grant access"
    # confirm button - STUDY_CATALOG is the same app-wide global the Study
    # Library itself lists from (global.R), not module-local state.
    approving_request_id <- reactiveVal(NULL)

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

    # Approve now opens a study picker rather than granting everything in one
    # click - admin chooses exactly which studies this request grants.
    observeEvent(input$approve_click, {
      approving_request_id(input$approve_click)
      showModal(modalDialog(
        title = "Grant catalog access",
        tags$p(class = "text-muted small", "Select which studies to grant access to."),
        checkboxGroupInput(session$ns("approve_study_ids"), NULL,
                           choices = setNames(STUDY_CATALOG$study_id, STUDY_CATALOG$display_name)),
        easyClose = TRUE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton(session$ns("confirm_approve"), "Grant access", class = "btn-success")
        )
      ))
    })

    observeEvent(input$confirm_approve, {
      rid <- approving_request_id()
      req(rid)
      result <- approve_request(db_path, rid, session$userData$user, input$approve_study_ids)
      showNotification(result$message, type = if (result$success) "message" else "error")
      if (result$success) {
        removeModal()
        approving_request_id(NULL)
      }
      refresh()
    })

    observeEvent(input$deny_click, {
      result <- deny_request(db_path, input$deny_click, session$userData$user)
      showNotification(result$message, type = if (result$success) "message" else "error")
      refresh()
    })
  })
}
