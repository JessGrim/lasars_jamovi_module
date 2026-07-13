# Shared helper functions for the lasars jamovi module
# These functions are called from lasarsmodel.b.R and origmodel.b.R.

# -----------------------------------------------------------------------------
# Input validation helpers
# -----------------------------------------------------------------------------

check_required_inputs <- function(options) {
  is.null(options$subject) || is.null(options$response) || is.null(options$subscale)
}

check_direction_requirements <- function(self,
                                         direction_required = TRUE,
                                         dir_var_name = 'direction',
                                         rev_label_name = 'rev_label') {

  if (!direction_required && is.null(self$options[[dir_var_name]]))
    return(TRUE)

  if (direction_required && is.null(self$options[[dir_var_name]])) {
    msg <- "This analysis requires a directional variable.\n" 
    msg <- paste0(msg,
                  "If none exists, use the Latent States and Response Style Analysis")
    notice <- jmvcore::Notice$new(options = self$options, name = 'dirMissing')
    notice$setContent(msg)
    self$results$insert(1, notice)
    return(FALSE)
  }

  if (!is.null(self$options[[dir_var_name]]) &&
      (is.null(self$options[[rev_label_name]]) ||
       is.na(self$options[[rev_label_name]]) ||
       self$options[[rev_label_name]] == "")) {
    msg <- paste0("Please enter the reverse-coding identifying label.\n",
                  "This is the string or character from the '",
                  self$options[[dir_var_name]], "' column\n",
                  "used to indicate that an item is reverse-coded.")
    notice <- jmvcore::Notice$new(options = self$options, name = 'negNeeded')
    notice$setContent(msg)
    self$results$insert(1, notice)
    return(FALSE)
  }

  if (!is.null(self$options[[dir_var_name]]) &&
      !is.null(self$options[[rev_label_name]]) &&
      self$options[[rev_label_name]] != "") {
    valid_levels <- unique(self$data[[self$options[[dir_var_name]]]])
    if (!(self$options[[rev_label_name]] %in% valid_levels)) {
      msg <- paste0("Could not find level '",
                    self$options[[rev_label_name]],
                    "' in ",
                    self$options[[dir_var_name]],
                    " column. Please enter a valid reverse-coding identifier.
                    Note: This variable is case sensitive.")
      notice <- jmvcore::Notice$new(options = self$options, name = 'invalidRevLabel')
      notice$setContent(msg)
      self$results$insert(1, notice)
      return(FALSE)
    }
  }

  TRUE
}

check_latent_states <- function(self) {
  scale_levels <- unique(self$data[[self$options$subscale]])
  if (!is.null(self$options$subscale) && length(scale_levels) <= 1) {
    msg <- paste0(
      "This analysis requires more than 1 latent state.\n",
      "The '", self$options$subscale, "' column only contains levels: ",
      paste(scale_levels, collapse = ", ")
    )
    notice <- jmvcore::Notice$new(options = self$options, name = 'tooFewLS')
    notice$setContent(msg)
    self$results$insert(1, notice)
    return(FALSE)
  }
  TRUE
}

check_response_options <- function(self) {
  resp_opts <- self$options$resp_opts

  if (is.null(resp_opts) || !is.numeric(resp_opts))
    return(TRUE)

  if (resp_opts <= 0) {
    msg <- "Number of response options must be greater than zero"
    notice <- jmvcore::Notice$new(options = self$options, name = 'invalidRespOpts')
    notice$setContent(msg)
    self$results$insert(1, notice)
    return(FALSE)
  }

  response_col <- self$options$response
  if (!is.null(response_col) && response_col != "") {
    response_values <- self$data[[response_col]]
    if (is.factor(response_values))
      response_values <- suppressWarnings(as.numeric(as.character(response_values)))

    if (any(!is.na(response_values) & response_values > resp_opts)) {
      msg <- paste0(
        "Values greater than ", resp_opts,
        " (current value in 'Number of response options' box) in ",
        response_col,
        " column"
      )
      notice <- jmvcore::Notice$new(options = self$options, name = 'outOfBoundsRespOpts')
      notice$setContent(msg)
      self$results$insert(1, notice)
      return(FALSE)
    }
  }

  TRUE
}

# -----------------------------------------------------------------------------
# Shared model utilities
# -----------------------------------------------------------------------------

make_thresholds <- function(K, centrePref_alpha, oddsPref_gamma, directionPref_omega) {
  lambda <- numeric(K - 1)
  G <- pnorm(oddsPref_gamma)
  lowers <- 1:floor(K / 2)
  lambda[lowers] <- exp(centrePref_alpha) * (-(K - 2) / 2) * G^(lowers - 1) - directionPref_omega
  uppers <- ceiling((K + 0.1) / 2):(K - 1)
  lambda[uppers] <- exp(centrePref_alpha) * ((K - 2) / 2) * G^(K - 1 - uppers) - directionPref_omega
  if ((K %% 2) == 0)
    lambda[K / 2] <- directionPref_omega
  lambda
}

build_mu_parameter_names <- function(data, subscale) {
  if (is.null(subscale) || is.null(data[[subscale]]))
    return(character())
  paste0("mu.", unique(data[[subscale]]))
}

build_lasars_response_style_names <- function(estDirectPref) {
  if (isTRUE(estDirectPref))
    c("centrePref", "oddsPref", "directionPref")
  else
    c("centrePref", "oddsPref")
}

build_subject_row_count <- function(data, subject, response = NULL) {
  if (!is.null(response) && !is.null(data[[response]]))
    data <- data[!is.na(data[[response]]), ]
  length(unique(data[[subject]]))
}

process_results <- function(results_list, data, subject, mu_indices = NULL,
                            upper_bound = 20, lower_bound = -20) {
  results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    if (is.null(mu_indices)) {
      c(
        round(results_list[[i]]$par, 3),
        convergence = results_list[[i]]$convergence
      )
    } else {
      c(
        round(results_list[[i]]$par[mu_indices], 3),
        convergence = results_list[[i]]$convergence
      )
    }
  }))
  results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)
  rownames(results_df) <- unique(data[[subject]])
  results_df$convergence <- ifelse(results_df$convergence == 0, "Success",
    ifelse(results_df$convergence == 1, "Error - Limit Reached",
      paste0("Error - optim code: ", results_df$convergence)))
  check_cols <- setdiff(names(results_df), "convergence")
  has_upper <- apply(results_df[check_cols], 1, function(x) any(x == upper_bound))
  has_lower <- apply(results_df[check_cols], 1, function(x) any(x == lower_bound))
  results_df$convergence <- ifelse(
    has_upper & has_lower, "Both bounds reached",
    ifelse(
      has_upper, "Upper Bound Reached",
      ifelse(
        has_lower, "Lower Bound Reached",
        results_df$convergence
      )
    )
  )
  results_df$subject <- rownames(results_df)
  results_df <- results_df[, c("subject", setdiff(names(results_df), "subject"))]
  results_df
}

ensure_table_rows <- function(table, row_count) {
  current_rows <- tryCatch({
    table$rowCount
  }, error = function(e) {
    0
  })
  if (is.na(current_rows) || current_rows < 0)
    current_rows <- 0
  if (row_count <= 0) {
    if (current_rows > 0)
      table$deleteRows()
    return()
  }
  if (current_rows > row_count) {
    table$deleteRows()
    current_rows <- 0
  }
  while (current_rows < row_count) {
    table$addRow(rowKey = current_rows + 1)
    current_rows <- current_rows + 1
  }
}

fill_subject_results_table <- function(results_table, results_df) {
  existing_cols <- character()
  if (!is.null(results_table$columns) && length(results_table$columns) > 0) {
    existing_cols <- vapply(results_table$columns, function(col) col$name, character(1))
  }
  for (add_col in setdiff(colnames(results_df), existing_cols)) {
    results_table$addColumn(name = add_col, format = "zto,3")
  }
  ensure_table_rows(results_table, nrow(results_df))
  for (i in seq_len(nrow(results_df))) {
    results_table$setRow(
      rowNo = i,
      values = as.list(results_df[i, , drop = FALSE])
    )
  }
  results_table$setVisible(visible = TRUE)
}

fill_descriptive_table <- function(desc_table, results_df, param_rows,
                                   note_key = "note", note_text = NULL) {
  if (length(param_rows) == 0)
    return()
  ensure_table_rows(desc_table, length(param_rows))
  for (i in seq_len(length(param_rows))) {
    desc_table$setRow(
      rowNo = i,
      values = list(
        par = param_rows[i],
        num = length(results_df[[param_rows[i]]]),
        mean = mean(results_df[[param_rows[i]]]),
        med = median(results_df[[param_rows[i]]]),
        sd = sd(results_df[[param_rows[i]]]),
        se = sd(results_df[[param_rows[i]]], na.rm = TRUE) / sqrt(sum(!is.na(results_df[[param_rows[i]]])))
      )
    )
  }
  if (!is.null(note_text))
    desc_table$setNote(key = note_key, note = note_text)
}

calculate_font_size <- function(x) {
  if (x == 0) {
    3
  } else {
    3 + 17 * (abs(x) - 0) / (1 - 0)
  }
}

format_correlation <- function(cor_value) {
  cor_char <- as.character(cor_value)
  if (substr(cor_char, 1, 1) == '-') {
    if (cor_char != "-1")
      cor_char <- paste0("-", substring(cor_char, first = 3))
  } else {
    if (cor_char != "1" && cor_char != "0")
      cor_char <- substring(cor_char, first = 2)
  }
  cor_char
}

build_pairs_plot <- function(plotData) {
  if (is.null(plotData) || ncol(plotData) < 2)
    return(NULL)
  just_scale_mus <- plotData[, 2:ncol(plotData), drop = FALSE]
  p_model <- GGally::ggpairs(plotData,
                             columns = 2:ncol(plotData),
                             upper = 'blank',
                             lower = list(continuous = GGally::wrap("points", size = 0.6, alpha = 1)),
                             axisLabels = "none",
                             columnLabels = NULL)
  for (i in seq_len(ncol(just_scale_mus))) {
    this_scale <- colnames(just_scale_mus)[i]
    p_model[i, i] <- GGally::ggally_text(this_scale,
                                        color = "black",
                                        size = 3.5)
  }
  if (ncol(just_scale_mus) > 1) {
    for (i in seq_len(ncol(just_scale_mus) - 1)) {
      for (j in (i + 1):ncol(just_scale_mus)) {
        this_cor <- round(cor(just_scale_mus[i], just_scale_mus[j], method = 'spearman'), 2)
        font_size <- as.numeric(calculate_font_size(this_cor))
        this_cor_char <- format_correlation(this_cor)
        p_model[i, j] <- GGally::ggally_text(this_cor_char,
                                             color = "black",
                                             size = font_size)
      }
    }
  }
  p_model + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = "white"),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA))
}
