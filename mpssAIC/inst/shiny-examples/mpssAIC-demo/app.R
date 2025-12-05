# inst/shiny-examples/mpssAIC-demo/app.R

library(shiny)
library(mpssAIC) # your package with mpssAIC:::multi_path_AIC_pipeline(), confusion_metrics(), etc.

## =========================
## UI
## =========================
ui <- fluidPage(
  titlePanel("mpssAIC – Multi-Path AIC Pipeline Demo"),
  tabsetPanel(
    id = "main_tabs",

    # -----------------------
    # Tab 1: Linear (Gaussian)
    # -----------------------
    tabPanel(
      "Linear Regression (Gaussian)",
      br(),
      h3("Plausible Models (Gaussian)"),
      tableOutput("linear_plausible"),
      br(),
      h3("Stability Scores"),
      tableOutput("linear_stability"),
      br(),
      h3("Confusion Matrix (Best Plausible Model)"),
      tableOutput("linear_cm"),
      br(),
      h3("Classification Metrics"),
      tableOutput("linear_metrics"),
      br(),
      h3("Confusion Matrix Plot"),
      plotOutput("linear_cm_plot", height = "300px")
    ),

    # -----------------------
    # Tab 2: Logistic (Binomial)
    # -----------------------
    tabPanel(
      "Logistic Regression (Binomial)",
      br(),
      h3("Plausible Models (Binomial)"),
      tableOutput("logistic_plausible"),
      br(),
      h3("Stability Scores"),
      tableOutput("logistic_stability"),
      br(),
      h3("Confusion Matrix (Best Plausible Model)"),
      tableOutput("logistic_cm"),
      br(),
      h3("Classification Metrics"),
      tableOutput("logistic_metrics"),
      br(),
      h3("Confusion Matrix Plot"),
      plotOutput("logistic_cm_plot", height = "300px")
    )
  )
)

## =========================
## SERVER
## =========================
server <- function(input, output, session) {
  ## --------------------------------------------------------------
  ## 1. Simulate synthetic data for both examples
  ## --------------------------------------------------------------
  set.seed(123)

  # Linear regression (Gaussian)
  n_lin <- 300
  p_lin <- 6
  X_linear <- matrix(rnorm(n_lin * p_lin), nrow = n_lin, ncol = p_lin)
  colnames(X_linear) <- paste0("x", 1:p_lin)

  beta_lin <- c(1.0, -1.5, 0.75, 0.5, 0, 0)
  y_linear <- as.numeric(X_linear %*% beta_lin + rnorm(n_lin, sd = 1))

  # Logistic regression (Binomial)
  n_log <- 300
  p_log <- 6
  X_logistic <- matrix(rnorm(n_log * p_log), nrow = n_log, ncol = p_log)
  colnames(X_logistic) <- paste0("x", 1:p_log)

  beta_log <- c(1.0, -1.5, 0.75, 0, 0, 0)
  eta <- X_logistic %*% beta_log
  prob <- plogis(eta)
  y_logistic <- rbinom(n_log, size = 1, prob = prob)

  ## --------------------------------------------------------------
  ## 2. Run multi-path AIC pipeline for Gaussian and Binomial
  ## --------------------------------------------------------------

  # Gaussian
  result_linear <- mpssAIC:::multi_path_AIC_pipeline(
    X             = X_linear,
    y             = y_linear,
    model_type    = "gaussian",
    K             = NULL,
    eps           = 1e-6,
    delta         = 2,
    L             = 50,
    B             = 20, # increase for more stable results
    resample_type = "bootstrap",
    m             = NULL,
    Delta         = 2,
    tau           = 0.6
  )

  # Binomial
  result_logistic <- mpssAIC:::multi_path_AIC_pipeline(
    X             = X_logistic,
    y             = y_logistic,
    model_type    = "binomial",
    K             = NULL,
    eps           = 1e-6,
    delta         = 2,
    L             = 50,
    B             = 20,
    resample_type = "bootstrap",
    m             = NULL,
    Delta         = 2,
    tau           = 0.6
  )

  ## --------------------------------------------------------------
  ## 3. Linear outputs: plausible models & stability
  ## --------------------------------------------------------------

  output$linear_plausible <- renderTable({
    if (is.null(result_linear$plausible) ||
      nrow(result_linear$plausible) == 0) {
      return(NULL)
    }
    tab <- result_linear$plausible[
      , c("AIC", "AIC_diff", "size", "vars", "avg_stability")
    ]
    rownames(tab) <- NULL
    tab
  })

  output$linear_stability <- renderTable({
    if (is.null(result_linear$stability_pi)) {
      return(NULL)
    }
    data.frame(
      Variable  = names(result_linear$stability_pi),
      Stability = round(as.numeric(result_linear$stability_pi), 3),
      row.names = NULL
    )
  })

  ## --------------------------------------------------------------
  ## 4. Linear: confusion matrix + metrics (best plausible model)
  ## --------------------------------------------------------------

  cm_results_linear <- reactive({
    if (is.null(result_linear$plausible) ||
      nrow(result_linear$plausible) == 0) {
      return(NULL)
    }

    # Best plausible model (first row)
    best_key <- result_linear$plausible$key[1]
    best_idx <- as.integer(strsplit(best_key, ",")[[1]])

    # Subset X and build data frame
    X_best <- X_linear[, best_idx, drop = FALSE]
    df_linear <- data.frame(y = y_linear, X_best)

    # Fit linear model
    fit_best <- lm(y ~ ., data = df_linear)

    # Binary outcome: high vs low (median split)
    cutoff_y <- median(y_linear, na.rm = TRUE)
    y_true <- as.integer(y_linear >= cutoff_y)

    # Predicted scores from linear model
    y_hat <- predict(fit_best, type = "response")

    # Rescale predictions to [0,1] to act like probabilities
    rng <- range(y_hat, na.rm = TRUE)
    if (diff(rng) == 0) {
      p_hat <- rep(0.5, length(y_hat))
    } else {
      p_hat <- (y_hat - rng[1]) / (rng[2] - rng[1])
    }

    cutoff <- 0.5
    y_pred <- as.integer(p_hat >= cutoff)

    # Confusion matrix for display/plot
    tab <- table(Actual = y_true, Predicted = y_pred)

    # Metrics using package confusion_metrics(y_true, p_hat, cutoff)
    metrics_df <- confusion_metrics(
      y_true = y_true,
      p_hat  = p_hat,
      cutoff = cutoff
    )

    metrics_vec <- as.numeric(metrics_df[1, ])
    names(metrics_vec) <- names(metrics_df)

    list(
      confusion_matrix = tab,
      metrics          = metrics_vec
    )
  })

  output$linear_cm <- renderTable(
    {
      res <- cm_results_linear()
      if (is.null(res)) {
        return(NULL)
      }
      res$confusion_matrix
    },
    rownames = TRUE
  )

  output$linear_metrics <- renderTable({
    res <- cm_results_linear()
    if (is.null(res)) {
      return(NULL)
    }
    data.frame(
      Metric = names(res$metrics),
      Value = round(as.numeric(res$metrics), 4),
      row.names = NULL
    )
  })

  output$linear_cm_plot <- renderPlot({
    res <- cm_results_linear()
    if (is.null(res)) {
      return(NULL)
    }
    cm <- res$confusion_matrix

    par(mar = c(5, 5, 4, 2))

    plot(1, 1,
      type = "n", xlim = c(0, 2), ylim = c(0, 2),
      xlab = "Predicted", ylab = "Actual",
      main = "Confusion Matrix (Linear – Best Plausible Model)",
      xaxt = "n", yaxt = "n", bty = "n"
    )

    rect(c(0, 1), c(0, 0), c(1, 2), c(1, 1),
      col = c("lightgreen", "lightcoral")
    )
    rect(c(0, 1), c(1, 1), c(1, 2), c(2, 2),
      col = c("lightcoral", "lightgreen")
    )

    text(0.5, 0.5, paste0("TN\n", cm[1, 1]), cex = 1.5, font = 2)
    text(1.5, 0.5, paste0("FP\n", cm[1, 2]), cex = 1.5, font = 2)
    text(0.5, 1.5, paste0("FN\n", cm[2, 1]), cex = 1.5, font = 2)
    text(1.5, 1.5, paste0("TP\n", cm[2, 2]), cex = 1.5, font = 2)

    axis(1, at = c(0.5, 1.5), labels = c("0", "1"))
    axis(2, at = c(0.5, 1.5), labels = c("0", "1"))
  })

  ## --------------------------------------------------------------
  ## 5. Logistic outputs: plausible models & stability
  ## --------------------------------------------------------------

  output$logistic_plausible <- renderTable({
    if (is.null(result_logistic$plausible) ||
      nrow(result_logistic$plausible) == 0) {
      return(NULL)
    }
    tab <- result_logistic$plausible[
      , c("AIC", "AIC_diff", "size", "vars", "avg_stability")
    ]
    rownames(tab) <- NULL
    tab
  })

  output$logistic_stability <- renderTable({
    if (is.null(result_logistic$stability_pi)) {
      return(NULL)
    }
    data.frame(
      Variable  = names(result_logistic$stability_pi),
      Stability = round(as.numeric(result_logistic$stability_pi), 3),
      row.names = NULL
    )
  })

  ## --------------------------------------------------------------
  ## 6. Logistic: confusion matrix + metrics (best plausible model)
  ## --------------------------------------------------------------

  cm_results <- reactive({
    if (is.null(result_logistic$plausible) ||
      nrow(result_logistic$plausible) == 0) {
      return(NULL)
    }

    best_key <- result_logistic$plausible$key[1]
    best_idx <- as.integer(strsplit(best_key, ",")[[1]])

    X_best <- X_logistic[, best_idx, drop = FALSE]
    df_logistic <- data.frame(y = y_logistic, X_best)

    # Fit best logistic model
    fit_best <- glm(y ~ ., data = df_logistic, family = binomial())

    cutoff <- 0.5
    p_hat <- predict(fit_best, type = "response")
    y_true <- y_logistic
    y_pred <- as.integer(p_hat >= cutoff)

    # Confusion matrix for display/plot
    tab <- table(Actual = y_true, Predicted = y_pred)

    # Metrics using package confusion_metrics(y_true, p_hat, cutoff)
    metrics_df <- confusion_metrics(
      y_true = y_true,
      p_hat  = p_hat,
      cutoff = cutoff
    )

    metrics_vec <- as.numeric(metrics_df[1, ])
    names(metrics_vec) <- names(metrics_df)

    list(
      confusion_matrix = tab,
      metrics          = metrics_vec
    )
  })

  output$logistic_cm <- renderTable(
    {
      res <- cm_results()
      if (is.null(res)) {
        return(NULL)
      }
      res$confusion_matrix
    },
    rownames = TRUE
  )

  output$logistic_metrics <- renderTable({
    res <- cm_results()
    if (is.null(res)) {
      return(NULL)
    }
    data.frame(
      Metric = names(res$metrics),
      Value = round(as.numeric(res$metrics), 4),
      row.names = NULL
    )
  })

  output$logistic_cm_plot <- renderPlot({
    res <- cm_results()
    if (is.null(res)) {
      return(NULL)
    }
    cm <- res$confusion_matrix

    par(mar = c(5, 5, 4, 2))

    plot(1, 1,
      type = "n", xlim = c(0, 2), ylim = c(0, 2),
      xlab = "Predicted", ylab = "Actual",
      main = "Confusion Matrix (Logistic – Best Plausible Model)",
      xaxt = "n", yaxt = "n", bty = "n"
    )

    rect(c(0, 1), c(0, 0), c(1, 2), c(1, 1),
      col = c("lightgreen", "lightcoral")
    )
    rect(c(0, 1), c(1, 1), c(1, 2), c(2, 2),
      col = c("lightcoral", "lightgreen")
    )

    text(0.5, 0.5, paste0("TN\n", cm[1, 1]), cex = 1.5, font = 2)
    text(1.5, 0.5, paste0("FP\n", cm[1, 2]), cex = 1.5, font = 2)
    text(0.5, 1.5, paste0("FN\n", cm[2, 1]), cex = 1.5, font = 2)
    text(1.5, 1.5, paste0("TP\n", cm[2, 2]), cex = 1.5, font = 2)

    axis(1, at = c(0.5, 1.5), labels = c("0", "1"))
    axis(2, at = c(0.5, 1.5), labels = c("0", "1"))
  })
}

## =========================
## RUN APP
## =========================
shinyApp(ui, server)
