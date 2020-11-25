tscausal.tsmodel.distribution <- function(object, actual, fitted = NULL, alpha = 0.05, ...) {
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
  summary_table <- data.table(type = c("Average", "Cumulative"),
    actual = c(mean(post_actual), sum(post_actual)),
    prediction = c(mean(point_pred), sum(point_pred)),
    prediction_lower = c(quantile(rowMeans(object), prob_lower), quantile(rowSums(object), prob_lower)),
    prediction_upper = c(quantile(rowMeans(object), prob_upper), quantile(rowSums(object), prob_upper)),
    prediction_sd = c(sd(rowMeans(object)),sd(rowSums(object))),
    absolute_effect = c(mean(post_actual) - mean(point_pred), sum(post_actual) - sum(point_pred)),
    absolute_effect_lower  = c(quantile(rowMeans(y_dist - object), prob_lower), quantile(rowSums(y_dist - object), prob_lower)),
    absolute_effect_upper = c(quantile(rowMeans(y_dist - object), prob_upper), quantile(rowSums(y_dist - object), prob_upper)),
    absolute_effect_sd = c(sd(rowMeans(y_dist - object)), sd(rowSums(y_dist - object))))
  summary_table[,relative_effect := absolute_effect/prediction]
  summary_table[,relative_effect_lower := absolute_effect_lower/prediction]
  summary_table[,relative_effect_upper := absolute_effect_upper/prediction]
  summary_table[,relative_effect_sd := absolute_effect_sd/prediction]
  # Add interval coverage, defined by alpha
  summary_table[,alpha := rep(alpha,2)]
  # Add one-sided tail-area probability of overall impact, p
  y_samples_post_sum <- rowSums(object)
  y_post_sum <- sum(post_actual)
  p <- min(sum(c(y_samples_post_sum, y_post_sum) >= y_post_sum),
           sum(c(y_samples_post_sum, y_post_sum) <= y_post_sum))/(length(y_samples_post_sum) + 1)
  if (p > 1 | p < 0) stop("\ntscausal-->error: estimated probability outside admissible range [0,1]...check.")
  summary_table[,p :=  p]
  out = list(summary_table = summary_table, predicted = object, post_actual = post_actual, pre_actual = pre_actual,
             fitted = fitted, intervention_date = intervention_date, alpha = alpha)
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
  colnames(f_summary) <- c("Average", "Cumulative")
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
    stmt <- NULL
    # Summarize averages
    stmt <- paste0(stmt, "\n\nDuring the post-intervention period, the response ",
                   "variable had an average value of approx. ", actual[1],
                   ". ", if (sig) "By contrast, in " else "In ",
                   "the absence of an intervention, ",
                   "we would have expected an average response of ",
                   prediction[1], ". The ", ci_coverage, " interval of this ",
                   "counterfactual prediction is [", prediction_lower[1],
                   ", ", prediction_upper[1], "]. Subtracting this ",
                   "prediction from the observed response yields an estimate ",
                   "of the causal effect the intervention had on the response ",
                   "variable. This effect is ", absolute_effect[1], " with a ",
                   ci_coverage, " interval of [", absolute_effect_lower[1], ", ",
                   absolute_effect_upper[1], "]. For a discussion of ",
                   "the significance of this effect, see below.")
    # Summarize sums
    stmt <- paste0(stmt, "\n\nSumming up the individual data points during ",
                   "the post-intervention period (which can only sometimes be ",
                   "meaningfully interpreted), the response variable had an ",
                   "overall value of ", actual[2], ". ",
                   if (sig) "By contrast, had " else "Had ",
                   "the intervention not taken place, we would have expected ",
                   "a sum of ", prediction[2], ". The ", ci_coverage, " interval of ",
                   "this prediction is [", prediction_lower[2], ", ", prediction_upper[2],
                   "].")
    # Summarize relative numbers (in which case row [1] = row [2])
    stmt <- paste0(stmt, "\n\nThe above results are given in terms of ",
                   "absolute numbers. In relative terms, the response variable ",
                   "showed ", if (pos) "an increase of " else "a decrease of",
                   relative_effect[1], ". The ", ci_coverage, " interval of this ",
                   "percentage is [", relative_effect_lower[1], ", ",
                   relative_effect_upper[1], "].")
    # Comment on significance
    if (sig && pos) {
      stmt <- paste0(stmt, "\n\nThis means that the positive effect observed ",
                     "during the intervention period is statistically ",
                     "significant and unlikely to be due to random ",
                     "fluctuations. ",
                     "It should be noted, however, that the question of whether ",
                     "this increase also bears substantive significance can ",
                     "only be answered by comparing the absolute effect (",
                     absolute_effect[1], ") to the original goal of ",
                     "the underlying intervention.")
    } else if (sig && !pos) {
      stmt <- paste0(stmt, "\n\nThis means that the negative effect observed ",
                     "during the intervention period is statistically ",
                     "significant. If the experimenter had expected a positive ",
                     "effect, it is recommended to ",
                     "double-check whether anomalies in the control variables ",
                     "may have caused an overly optimistic expectation of ",
                     "what should have happened in the response variable in the ",
                     "absence of the intervention.")
    } else if (!sig && pos) {
      stmt <- paste0(stmt, "\n\nThis means that, although the intervention ",
                     "appears to have caused a positive effect, this effect ",
                     "is not statistically significant when considering the ",
                     "entire post-intervention period as a whole. Individual ",
                     "days or shorter stretches within the intervention period ",
                     "may of course still have had a significant effect, as ",
                     "indicated whenever the lower limit of the impact ",
                     "time series (lower plot) was above zero.")
    } else if (!sig && !pos) {
      stmt <- paste0(stmt, "\n\nThis means that, although it may look as ",
                     "though the intervention has exerted a negative effect ",
                     "on the response variable when considering the ",
                     "intervention period as a whole, this effect is not ",
                     "statistically significant, and so cannot be ",
                     "meaningfully interpreted.")
    }
    if (!sig) {
      stmt <- paste0(stmt, " The apparent effect could be the result of ",
                     "random fluctuations that are unrelated to the ",
                     "intervention. This is often the case when the ",
                     "intervention period is very long and includes much ",
                     "of the time when the effect has already worn off. ",
                     "It can also be the case when the intervention period ",
                     "is too short to distinguish the signal from the noise. ",
                     "Finally, failing to find a significant effect can ",
                     "happen when there are not enough control variables or ",
                     "when these variables do not correlate well with ",
                     "the response variable during the learning period.")
    }
    if (p < summary_table$alpha[1]) {
      stmt <- paste0(stmt, "\n\nThe probability of obtaining this effect by ",
                     "chance is very small (Bayesian one-sided tail-area ",
                     "probability p = ", round(p, 3), "). This means the causal ",
                     "effect can be considered statistically significant.")
    } else {
      stmt <- paste0(stmt, "\n\nThe probability of obtaining this ",
                     "effect by chance is p = ", round(p, 3), ". This ",
                     "means the effect may be spurious and would generally ",
                     "not be considered statistically significant.")
    }
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
  # Initialize statement
  stmt <- NULL
  stmt <- cat(stmt, "  \n During the post-intervention period, the response ",
              "variable had an average value of approx. ", actual[1],
              ". ", if (sig) "By contrast, in " else "In ",
              "the absence of an intervention, ",
              "we would have expected an average response of ",
              prediction[1], ". The ", ci_coverage, " interval of this ",
              "counterfactual prediction is [", prediction_lower[1],
              ", ", prediction_upper[1], "]. Subtracting this ",
              "prediction from the observed response yields an estimate ",
              "of the causal effect the intervention had on the response ",
              "variable. This effect is ", absolute_effect[1], " with a ",
              ci_coverage, " interval of [", absolute_effect_lower[1], ", ",
              absolute_effect_upper[1], "]. For a discussion of ",
              "the significance of this effect, see below.")
  # Summarize sums
  stmt <- cat(stmt, "  \n Summing up the individual data points during ",
              "the post-intervention period (which can only sometimes be ",
              "meaningfully interpreted), the response variable had an ",
              "overall value of ", actual[2], ". ",
              if (sig) "By contrast, had " else "Had ",
              "the intervention not taken place, we would have expected ",
              "a sum of ", prediction[2], ". The ", ci_coverage, " interval of ",
              "this prediction is [", prediction_lower[2], ", ", prediction_upper[2],
              "].")
  # Summarize relative numbers (in which case row [1] = row [2])
  stmt <- cat(stmt,"  \n The above results are given in terms of ",
              "absolute numbers. In relative terms, the response variable ",
              "showed ", if (pos) "an increase of " else "a decrease of",
              relative_effect[1], ". The ", ci_coverage, " interval of this ",
              "percentage is [", relative_effect_lower[1], ", ",
              relative_effect_upper[1], "].")
  # Comment on significance
  if (sig && pos) {
    stmt <- cat(stmt, "  \n This means that the positive effect observed ",
                "during the intervention period is statistically ",
                "significant and unlikely to be due to random ",
                "fluctuations. ",
                "It should be noted, however, that the question of whether ",
                "this increase also bears substantive significance can ",
                "only be answered by comparing the absolute effect (",
                absolute_effect[1], ") to the original goal of ",
                "the underlying intervention.")
  } else if (sig && !pos) {
    stmt <- cat(stmt, "  \n This means that the negative effect observed ",
                "during the intervention period is statistically ",
                "significant. If the experimenter had expected a positive ",
                "effect, it is recommended to ",
                "double-check whether anomalies in the control variables ",
                "may have caused an overly optimistic expectation of ",
                "what should have happened in the response variable in the ",
                "absence of the intervention.")
  } else if (!sig && pos) {
    stmt <- cat(stmt, "  \n This means that, although the intervention ",
                "appears to have caused a positive effect, this effect ",
                "is not statistically significant when considering the ",
                "entire post-intervention period as a whole. Individual ",
                "days or shorter stretches within the intervention period ",
                "may of course still have had a significant effect, as ",
                "indicated whenever the lower limit of the impact ",
                "time series (lower plot) was above zero.")
  } else if (!sig && !pos) {
    stmt <- cat(stmt, "  \n This means that, although it may look as ",
                "though the intervention has exerted a negative effect ",
                "on the response variable when considering the ",
                "intervention period as a whole, this effect is not ",
                "statistically significant, and so cannot be ",
                "meaningfully interpreted.")
  }
  if (!sig) {
    stmt <- cat(stmt, " The apparent effect could be the result of ",
                "random fluctuations that are unrelated to the ",
                "intervention. This is often the case when the ",
                "intervention period is very long and includes much ",
                "of the time when the effect has already worn off. ",
                "It can also be the case when the intervention period ",
                "is too short to distinguish the signal from the noise. ",
                "Finally, failing to find a significant effect can ",
                "happen when there are not enough control variables or ",
                "when these variables do not correlate well with ",
                "the response variable during the learning period.")
  }
  if (p < summary_table$alpha[1]) {
    stmt <- cat(stmt, "  \n The probability of obtaining this effect by ",
                "chance is very small (Bayesian one-sided tail-area ",
                "probability p = ", round(p, 3), "). This means the causal ",
                "effect can be considered statistically significant.")
  } else {
    stmt <- cat(stmt, "  \n The probability of obtaining this ",
                "effect by chance is p = ", round(p, 3), ". This ",
                "means the effect may be spurious and would generally ",
                "not be considered statistically significant.")
  }
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
  colnames(f_summary) <- c("Average", "Cumulative")
  ci_label <- paste0(round((1 - alpha) * 100), "% CI")
  row.names(f_summary) <- c("Actual", "Prediction (s.d.)", ci_label,
                            " ",
                            "Absolute effect (s.d.)", paste(ci_label, ""),
                            "  ",
                            "Relative effect (s.d.)", paste(ci_label, " "))
  return(f_summary)
}

