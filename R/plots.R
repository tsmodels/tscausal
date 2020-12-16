plot.tscausal <- function(x, y = NULL, median_color = "black", median_type = 1, median_width = 2, 
                          gradient_color = "orange", interval_color = "red", interval_type = 2, 
                          interval_width = 2, ylim = NULL, ylab = "", n_original = NULL, ...)
{
  p <- process.tscausal(x, n_original = n_original)
  if (is.null(n_original)) n_original <- NROW(x$pre_actual)
  par(bg = "white", mar = c(1,2,0.5,3))
  if (x$include_cumulative) {
    layout(mat = matrix(c(1,2,3), nrow = 3))
    plot(p$prediction, interval_quantiles = c(p$alpha/2, 1 - p$alpha/2), interval_color = interval_color, gradient_color = gradient_color, x_axes = FALSE, 
         median_color = median_color, median_type = median_type, median_width = median_width, interval_type = interval_type, 
         interval_width = interval_width, ylim = ylim, ylab = ylab)
    abline(v = as.Date(p$intervention_date), col = "green", lty = 2, lwd = 2)
    mtext("Original", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    legend("topleft", c("Actual","Counterfactual"), col = c("red","black"), lty = 1, bty = "n", lwd = c(1,2))
    p$point_distribution = na.fill(p$point_distribution, fill = 0)
    class(p$point_distribution) <- "tsmodel.distribution"
    plot(p$point_distribution, interval_quantiles = c(p$alpha/2, 1 - p$alpha/2),interval_color = interval_color, gradient_color = gradient_color, x_axes = FALSE, 
         median_color = median_color, median_type = median_type, median_width = median_width, interval_type = interval_type, 
         interval_width = interval_width, ylim = ylim, ylab = ylab)
    abline(v = as.Date(p$intervention_date), col = "green", lty = 2, lwd = 2)
    mtext("Pointwise", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    par(mar = c(3,2,0.5,3))
    p$cumulative_effect_dist = na.fill(p$cumulative_effect_dist, fill = 0)
    class(p$cumulative_effect_dist) <- "tsmodel.distribution"
    plot(p$cumulative_effect_dist, interval_quantiles = c(p$alpha/2, 1 - p$alpha/2), interval_color = interval_color, gradient_color = gradient_color, x_axes = TRUE, 
         median_color = median_color, median_type = median_type, median_width = median_width, interval_type = interval_type, 
         interval_width = interval_width, ylim = ylim, ylab = ylab)
    abline(v = as.Date(p$intervention_date), col = "green", lty = 2, lwd = 2)
    mtext("Cumulative", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    layout(mat = matrix(c(1,2), nrow = 2))
    plot(p$prediction, interval_quantiles = c(p$alpha/2, 1 - p$alpha/2), interval_color = interval_color, gradient_color = gradient_color, x_axes = FALSE, 
         median_color = median_color, median_type = median_type, median_width = median_width, interval_type = interval_type, 
         interval_width = interval_width, ylim = ylim, ylab = ylab)
    abline(v = as.Date(p$intervention_date), col = "green", lty = 2, lwd = 2)
    mtext("Original", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    legend("topleft", c("Actual","Counterfactual"), col = c("red","black"), lty = 1, bty = "n", lwd = c(1,2))
    p$point_distribution = na.fill(p$point_distribution, fill = 0)
    class(p$point_distribution) <- "tsmodel.distribution"
    par(mar = c(3,2,0.5,3))
    plot(p$point_distribution, interval_quantiles = c(p$alpha/2, 1 - p$alpha/2), interval_color = interval_color, gradient_color = gradient_color, x_axes = TRUE, 
         median_color = median_color, median_type = median_type, median_width = median_width, interval_type = interval_type, 
         interval_width = interval_width, ylim = ylim, ylab = ylab)
    abline(v = as.Date(p$intervention_date), col = "green", lty = 2, lwd = 2)
    mtext("Pointwise", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
  }
}



process.tscausal <- function(object, n_original = NULL){
  # Check input
  alpha <- object$alpha
  prob_lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob_upper <- 1 - alpha / 2  # e.g., 0.975 when alpha = 0.05
  if (!is.null(object$fitted)) {
    if (is(object$fitted,"tsmodel.distribution")) {
      if (!is.null(n_original)) {
        N <- ncol(object$fitted)
        n <- pmin(pmax(1, n_original), N)
        fx <- object$fitted[,(N - n):N]
      } else{
        fx <- object$fitted
      }
      fit_dist <- fx
    } else{
      if (!is.null(n_original)) {
        fx <- tail(object$fitted, pmax(1, n_original))
      } else {
        fx <- object$fitted
      }
      fit_dist <- matrix(as.numeric(fx), ncol = NROW(fx), nrow = NROW(object$predicted), byrow = TRUE)
    }
  } else {
    if (is.null(n_original)) n_original <- NROW(object$pre_actual)
    fit_dist <- matrix(NA, ncol = n_original, nrow = NROW(object$predicted), byrow = TRUE)
  }
  colnames(fit_dist) <- as.character(tail(index(object$pre_actual), NCOL(fit_dist)))
  response <- rbind(tail(object$pre_actual, NCOL(fit_dist)), object$post_actual)
  total_pred_dist <- cbind(fit_dist, object$predicted)
  class(total_pred_dist) <- "tsmodel.distribution"
  actual_dist <- cbind(matrix(as.numeric(object$pre_actual[colnames(fit_dist)]), ncol = ncol(fit_dist), nrow = nrow(fit_dist), byrow = TRUE),
                      matrix(as.numeric(object$post_actual), ncol = NROW(object$post_actual), nrow = nrow(fit_dist),byrow = TRUE))
  colnames(actual_dist) <- c(colnames(fit_dist), as.character(index(object$post_actual)))
  class(actual_dist) <- "tsmodel.distribution"
  intN <- which(colnames(actual_dist) == object$intervention_date)
  N <- ncol(total_pred_dist)
  point_dist <- actual_dist - total_pred_dist
  class(point_dist) <- "tsmodel.distribution"
  pdist <- point_dist
  pdist[,1:intN] <- 0
  cumulative_effect <- t(apply(pdist, 1, cumsum))
  class(cumulative_effect) <- "tsmodel.distribution"
  pout <- list(distribution = object$predicted, original_series = zoo(rbind(object$pre_actual[colnames(fit_dist)], object$post_actual)), 
               h = ncol(object$predicted))
  class(pout) <- "tsmodel.predict"
  out <- list(prediction = pout, point_distribution = point_dist, cumulative_effect_dist = cumulative_effect, 
              intervention_date = object$intervention_date, alpha = object$alpha)
  return(out)
}
