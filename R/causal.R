tscausal.tsmodel.distribution <- function(object, actual, fitted = NULL, alpha = 0.05, include_cumulative = TRUE, ...) {
  # Check input
  pindices <- as.Date(colnames(object))
  post_actual <- actual[pindices]
  ix <- which(index(actual) %in% pindices)
  pre_actual <- actual[1:(ix[1] - 1)]
  intervention_date <- tail(index(pre_actual),1)
  if (!is.null(fitted)) pre_actual <- pre_actual[index(fitted)]
  if (NROW(post_actual) != ncol(object)) {
    stop("\ntscausal->error: length of y not equal to number of columns of forecast distribution (object).")
  }
  # We will compare the matrix of predicted trajectories (e.g., 900 x 201)
  # with a matrix of replicated observations (e.g., 900 x 201)
  y_dist <- matrix(as.numeric(post_actual), nrow = nrow(object), ncol = ncol(object), byrow = TRUE)
  point_pred <- colMeans(object)
  # Define quantiles
  alpha <- alpha[1]
  if (alpha > 1 | alpha < 0) stop("\ntscausal: alpha (coverage) must be between 0 and 1")
  prob_lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob_upper <- 1 - alpha / 2  # e.g., 0.975 when alpha = 0.05
  # Compile summary statistics
  if (include_cumulative) {
    rep_n <- 2
    col_names <- c("Average", "Cumulative")
    summary_table <- data.table(type = col_names,
                                actual = c(mean(post_actual), sum(post_actual)),
                                prediction = c(mean(point_pred), sum(point_pred)),
                                prediction_lower = c(quantile(rowMeans(object), prob_lower), quantile(rowSums(object), prob_lower)),
                                prediction_upper = c(quantile(rowMeans(object), prob_upper), quantile(rowSums(object), prob_upper)),
                                prediction_sd = c(sd(rowMeans(object)),sd(rowSums(object))),
                                absolute_effect = c(mean(post_actual) - mean(point_pred), sum(post_actual) - sum(point_pred)),
                                absolute_effect_lower  = c(quantile(rowMeans(y_dist - object), prob_lower), quantile(rowSums(y_dist - object), prob_lower)),
                                absolute_effect_upper = c(quantile(rowMeans(y_dist - object), prob_upper), quantile(rowSums(y_dist - object), prob_upper)),
                                absolute_effect_sd = c(sd(rowMeans(y_dist - object)), sd(rowSums(y_dist - object))))
  } else {
    rep_n <- 1
    col_names <- "Average"
    summary_table <- data.table(type = col_names,
                                actual = mean(post_actual),
                                prediction = mean(point_pred),
                                prediction_lower = quantile(rowMeans(object), prob_lower),
                                prediction_upper = quantile(rowMeans(object), prob_upper),
                                prediction_sd = sd(rowMeans(object)),
                                absolute_effect = mean(post_actual) - mean(point_pred),
                                absolute_effect_lower  = quantile(rowMeans(y_dist - object), prob_lower),
                                absolute_effect_upper = quantile(rowMeans(y_dist - object), prob_upper),
                                absolute_effect_sd = sd(rowMeans(y_dist - object)))
  }
  summary_table[,relative_effect := absolute_effect/prediction]
  summary_table[,relative_effect_lower := absolute_effect_lower/prediction]
  summary_table[,relative_effect_upper := absolute_effect_upper/prediction]
  summary_table[,relative_effect_sd := absolute_effect_sd/prediction]
  # Add interval coverage, defined by alpha
  summary_table[,alpha := rep(alpha,rep_n)]
  # Add one-sided tail-area probability of overall impact, p
  y_samples_post_sum <- rowSums(object)
  y_post_sum <- sum(post_actual)
  p <- min(sum(c(y_samples_post_sum, y_post_sum) >= y_post_sum),
           sum(c(y_samples_post_sum, y_post_sum) <= y_post_sum))/(length(y_samples_post_sum) + 1)
  if (p > 1 | p < 0) stop("\ntscausal-->error: estimated probability outside admissible range [0,1]...check.")
  summary_table[,p :=  p]
  out = list(summary_table = summary_table, predicted = object, post_actual = post_actual, pre_actual = pre_actual,
             fitted = fitted, intervention_date = intervention_date, alpha = alpha, include_cumulative = include_cumulative)
  class(out) = "tscausal"
  return(out)
}

print.tscausal = function(x, digits = 4, ...)
{
  summary_table <- x$summary_table
  alpha <- x$alpha
  # Print title
  cat("Predictive inference {tscausal}\n")
  if (is.null(summary_table)) {
    cat("(Inference aborted)\n")
    return(invisible(NULL))
  }
  # Compile data frame with formatted numbers
  if (x$include_cumulative) {
    sep <- c("", "")
  } else {
    sep <- ""
  }
  f_summary <- data.table(
    actual = format_number(summary_table$actual, digits = digits),
    prediction = paste0(format_number(summary_table$prediction, digits = digits), " (", format_number(summary_table$prediction_sd, digits = digits), ")"),
    prediction_ci = format_confidence_intervals(summary_table$prediction_lower, summary_table$prediction_upper, digits = digits),
    Separator1 = sep,
    absolute_effect = paste0(format_number(summary_table$absolute_effect, digits = digits)," (", format_number(summary_table$absolute_effect_sd, digits = digits), ")"),
    absolute_effect_ci = format_confidence_intervals(summary_table$absolute_effect_lower, summary_table$absolute_effect_upper, digits = digits),
    Separator2 = sep,
    relative_effect = paste0(format_percent(summary_table$relative_effect, digits = digits), " (", format_percent(summary_table$relative_effect_sd, digits = digits), ")"),
    relative_effect_ci = format_percent_confidence_intervals(summary_table$relative_effect_lower, summary_table$relative_effect_upper, digits = digits))
  # Invert and format as table
  f_summary <- t(f_summary)
  if (x$include_cumulative) {
    colnames(f_summary) <- c("Average","Cumulative")
  } else {
    colnames(f_summary) <- "Average"
  }
  ci_label <- paste0(round((1 - alpha) * 100), "% CI")
  row.names(f_summary) <- c("Actual", "Prediction (s.d.)", ci_label,
                            " ",
                            "Absolute effect (s.d.)", paste(ci_label, ""),
                            "  ",
                            "Relative effect (s.d.)", paste(ci_label, " "))
  cat("\n")
  print.default(f_summary, print.gap = 3L, quote = FALSE)
  cat("\n")
  # Print overall tail-area probability
  p <- summary_table$p[1]
  cat(paste0("Predictive Distribution tail-area probability p:   ", round(p, 5), "\n"))
  cat(paste0("Predictive Distribution prob. of a causal effect:  ", round((1 - p) * 100, ifelse(p < 0.01, 5, ifelse(p < 0.05, 3, 0))), "%\n"))
  cat("\n")
}

# report statement summary output is based on CausalImpace Code
tsreport.tscausal = function(object, digits = 4, doc_template = NULL, type = c("screen", "pdf", "doc", "html"), output_dir = "/", args = list(name = "Causal Analysis", frequency = NULL, model = NULL), ...)
{
  type <- match.arg(type[1], c("screen", "pdf", "doc", "html"))
  if (type == "screen") {
    summary_table = object$summary_table
    actual <- prettify_number(summary_table$actual, round_digits = digits)
    letter <- identify_number_abbreviation(actual)
    prediction <- prettify_number(summary_table$prediction, letter, 2)
    prediction_lower <- prettify_number(summary_table$prediction_lower, letter, digits)
    prediction_upper <- prettify_number(summary_table$prediction_upper, letter, digits)
    absolute_effect <- prettify_number(summary_table$absolute_effect, letter, digits)
    absolute_effect_lower <- prettify_number(summary_table$absolute_effect_lower, letter, digits)
    absolute_effect_upper <- prettify_number(summary_table$absolute_effect_upper, letter, digits)
    relative_effect <- prettify_percentage(summary_table$relative_effect)
    relative_effect_lower <- prettify_percentage(summary_table$relative_effect_lower)
    relative_effect_upper <- prettify_percentage(summary_table$relative_effect_upper)
    # Evaluate significance and direction of the effect (increase or decrease)
    sig <- (!((summary_table$relative_effect_lower[1] < 0) && (summary_table$relative_effect_upper[1] > 0)))
    pos <- summary_table$relative_effect[1] > 0
    p <- summary_table$p[1]
    # Interval name
    ci_coverage <- paste0(round((1 - summary_table$alpha[1]) * 100), "%")
    # Initialize statement
    # Summarize averages
    stmt <- report_statement(actual, prediction, object, summary_table$alpha[1], p, sig, pos, ci_coverage, prediction_lower, prediction_upper, absolute_effect, 
                                 relative_effect, relative_effect_lower, relative_effect_upper, 
                                 absolute_effect_lower, absolute_effect_upper)
    return(stmt)
  } else {
    saveRDS(object, file = paste0(output_dir,"/causal_tmp.rds"))
    fname <- args$name
    if (is.null(fname)) fname <- paste0(fname,"_tscausal_report_",format(as.Date(Sys.Date()),"%Y%m%d")) else fname = paste0("tscausal_report_",format(as.Date(Sys.Date()),"%Y%m%d"))
    if (type == "pdf") {
      suppressWarnings(rmarkdown::render(file.path(find.package('tscausal'),'scripts/causal_report_pdf.Rmd'),
                                         output_dir = output_dir,
                                         output_file = paste0(fname,".pdf"),
                                         intermediates_dir = tempdir(),
                                         output_format = pdf_document2(toc = FALSE),
                                         params = list(dir = as.character(output_dir), name = args$name, frequency = args$frequency, model = args$model), quiet = TRUE))
    } else if (type == "doc") {
      file.copy(file.path(find.package('tscausal'),'scripts/causal_report_docx.Rmd'), file.path(paste0(output_dir,"/causal_report_docx.Rmd")))
      suppressWarnings(rmarkdown::render(file.path(paste0(output_dir,"/causal_report_docx.Rmd")),
                                         output_dir = output_dir,
                                         output_file = paste0(fname,".docx"),
                                         intermediates_dir = tempdir(),
                                         output_format = word_document2(toc = FALSE, reference_docx = doc_template, title = "Time Series Causal Analysis"),
                                         params = list(dir = as.character(output_dir), name = args$name, frequency = args$frequency, model = args$model), quiet = F))
      
    } else {
      suppressWarnings(rmarkdown::render(file.path(find.package('tscausal'),'scripts/causal_report_html.Rmd'),
                                         output_dir = output_dir,
                                         output_file = paste0(fname,".html"),
                                         intermediates_dir = tempdir(),
                                         output_format = html_document2(toc = FALSE),
                                         params = list(dir = as.character(output_dir), name = args$name, frequency = args$frequency, model = args$model), quiet = TRUE))
    }
  }
}

.latex.output = function(x, digits = 4, ...)
{
  summary_table <- x$summary_table
  actual <- prettify_number(summary_table$actual, round_digits = digits)
  letter <- identify_number_abbreviation(actual)
  prediction <- prettify_number(summary_table$prediction, letter, 2)
  prediction_lower <- prettify_number(summary_table$prediction_lower, letter, digits)
  prediction_upper <- prettify_number(summary_table$prediction_upper, letter, digits)
  absolute_effect <- prettify_number(summary_table$absolute_effect, letter, digits)
  absolute_effect_lower <- prettify_number(summary_table$absolute_effect_lower, letter, digits)
  absolute_effect_upper <- prettify_number(summary_table$absolute_effect_upper, letter, digits)
  relative_effect <- prettify_percentage(summary_table$relative_effect)
  relative_effect_lower <- prettify_percentage(summary_table$relative_effect_lower)
  relative_effect_upper <- prettify_percentage(summary_table$relative_effect_upper)
  # Evaluate significance and direction of the effect (increase or decrease)
  sig <- (!((summary_table$relative_effect_lower[1] < 0) && (summary_table$relative_effect_upper[1] > 0)))
  pos <- summary_table$relative_effect[1] > 0
  p <- summary_table$p[1]
  # Interval name
  ci_coverage <- paste0(round((1 - summary_table$alpha[1]) * 100), "%")
  stmt <- report_statement(actual, prediction, object = x, summary_table$alpha[1], p, sig, pos, ci_coverage, prediction_lower, prediction_upper, absolute_effect, 
                           relative_effect, relative_effect_lower, relative_effect_upper, 
                           absolute_effect_lower, absolute_effect_upper)
  stmt <- cat(stmt, " ")
  return(stmt)
}

.table.print.tscausal = function(x, digits = 4, ...)
{
  summary_table <- x$summary_table
  alpha <- x$alpha
  # Print title
  # Compile data frame with formatted numbers
  f_summary <- data.table(
    actual = format_number(summary_table$actual, digits = digits),
    prediction = paste0(format_number(summary_table$prediction, digits = digits), " (", format_number(summary_table$prediction_sd, digits = digits), ")"),
    prediction_ci = format_confidence_intervals(summary_table$prediction_lower, summary_table$prediction_upper, digits = digits),
    Separator1 = c("", ""),
    absolute_effect = paste0(format_number(summary_table$absolute_effect, digits = digits)," (", format_number(summary_table$absolute_effect_sd, digits = digits), ")"),
    absolute_effect_ci = format_confidence_intervals(summary_table$absolute_effect_lower, summary_table$absolute_effect_upper, digits = digits),
    Separator2 = c("", ""),
    relative_effect = paste0(format_percent(summary_table$relative_effect, digits = digits), " (", format_percent(summary_table$relative_effect_sd, digits = digits), ")"),
    relative_effect_ci = format_percent_confidence_intervals(summary_table$relative_effect_lower, summary_table$relative_effect_upper, digits = digits))
  # Invert and format as table
  f_summary <- t(f_summary)
  if (x$include_cumulative) {
    colnames(f_summary) <- c("Average", "Cumulative")
  } else {
    colnames(f_summary) <- "Average"
  }
  ci_label <- paste0(round((1 - alpha) * 100), "% CI")
  row.names(f_summary) <- c("Actual", "Prediction (s.d.)", ci_label,
                            " ",
                            "Absolute effect (s.d.)", paste(ci_label, ""),
                            "  ",
                            "Relative effect (s.d.)", paste(ci_label, " "))
  return(f_summary)
}