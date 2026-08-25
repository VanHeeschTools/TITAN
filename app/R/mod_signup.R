## Self-service sign-up module (email + password + confirm password, plus an
## optional catalog-access request).
##
## shinymanager has no built-in sign-up form, so this is injected into the
## login screen's "Sign up" pane (see .auth_tabs_close in app.R, which wraps
## this in a <div id="titan-tab-signup-pane"> toggled by .auth_tab_toggle_js).
## New accounts always get role = 'general'; checking "I need catalog access"
## immediately files a request via submit_catalog_request() (R/catalog_access.R
## — the same function/table the post-login "request access" flow uses), so
## it shows up in the admin's pending-requests queue (R/mod_admin_requests.R)
## right away instead of requiring a second trip through that flow.

mod_signup_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # autocomplete is patched onto these fields via JS (.auth_autocomplete_js
    # in app.R) - this Shiny version's textInput()/passwordInput() validate
    # `...` as empty (no pass-through HTML attrs), and tagAppendAttributes()
    # would land on the outer wrapper div, not the actual <input>.
    textInput(ns("email"), "Email:", width = "100%"),
    passwordInput(ns("password"), "Password (min 8 characters):", width = "100%"),
    passwordInput(ns("password_confirm"), "Confirm password:", width = "100%"),
    tags$hr(),
    checkboxInput(ns("want_catalog_access"),
                  "I need catalog access (browse the study library)",
                  value = FALSE),
    textAreaInput(ns("justification"), "Why do you need catalog access? (optional)",
                  rows = 2, width = "100%",
                  placeholder = "e.g. study/team, what you'll use it for"),
    tags$div(
      style = "text-align:center;",
      actionButton(ns("submit"), "Sign up", width = "100%", class = "btn-primary"),
      tags$br(), tags$br()
    ),
    tags$div(id = ns("result"))
  )
}

mod_signup_server <- function(id, db_path) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$submit, {
      email    <- trimws(input$email)
      password <- input$password
      confirm  <- input$password_confirm

      if (!grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", email, perl = TRUE)) {
        showNotification("Enter a valid email address.", type = "error")
        return(invisible(NULL))
      }
      if (nchar(password) < 8) {
        showNotification("Password must be at least 8 characters.", type = "error")
        return(invisible(NULL))
      }
      if (!identical(password, confirm)) {
        showNotification("Passwords do not match.", type = "error")
        return(invisible(NULL))
      }

      existing <- safe_db_read(db_path, "SELECT 1 FROM users WHERE email = ?", params = list(email))
      if (nrow(existing) > 0) {
        showNotification("An account with that email already exists.", type = "error")
        return(invisible(NULL))
      }

      result <- safe_db_write(
        db_path,
        "INSERT INTO users (email, password_hash, role) VALUES (?, ?, 'general')",
        params = list(email, hash_password(password))
      )

      if (isFALSE(result)) {
        showNotification("Could not create account — please try again.", type = "error")
        return(invisible(NULL))
      }

      account_msg <- "Account created — you can now log in."
      if (isTRUE(input$want_catalog_access)) {
        justification <- trimws(input$justification)
        req_result <- submit_catalog_request(
          db_path, email,
          justification = if (nzchar(justification)) justification else NA_character_
        )
        account_msg <- paste(account_msg, req_result$message)
      }

      showNotification(account_msg, type = "message", duration = 8)
      updateTextInput(session, "email", value = "")
      updateTextInput(session, "password", value = "")
      updateTextInput(session, "password_confirm", value = "")
      updateCheckboxInput(session, "want_catalog_access", value = FALSE)
      updateTextAreaInput(session, "justification", value = "")
    })
  })
}
