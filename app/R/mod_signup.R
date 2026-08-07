## Self-service sign-up module (email + password + confirm password).
##
## shinymanager has no built-in sign-up form, so this is injected into the
## login screen via secure_app(tags_bottom = ...) in app.R, hidden by default
## and toggled against the login fields with plain JS (see .signup_toggle_js
## in app.R) — no shinyjs dependency needed for a single show/hide toggle.
## New accounts always get role = 'general'.

mod_signup_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    id = "titan-signup-panel", style = "display:none; text-align:left;",
    tags$br(),
    textInput(ns("email"), "Email:", width = "100%"),
    passwordInput(ns("password"), "Password (min 8 characters):", width = "100%"),
    passwordInput(ns("password_confirm"), "Confirm password:", width = "100%"),
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

      if (!grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", email)) {
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

      showNotification("Account created — you can now log in.", type = "message", duration = 6)
      updateTextInput(session, "email", value = "")
      updateTextInput(session, "password", value = "")
      updateTextInput(session, "password_confirm", value = "")
    })
  })
}
