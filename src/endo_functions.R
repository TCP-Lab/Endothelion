

# Colored echo for on-screen log in Bash CLI
echo <- function(text, color = "white") {
  # Check the operating system
  os_type <- Sys.info()["sysname"]
  if (os_type == "Windows") {
    warning("Colored echo is just for Bash!")
    cat(text)
  } else if (os_type == "Linux") {
    if      (color == "red")     {col <- "\\e[1;31m"}
    else if (color == "green")   {col <- "\\e[1;32m"}
    else if (color == "yellow")  {col <- "\\e[1;33m"}
    else if (color == "blue")    {col <- "\\e[1;34m"}
    else if (color == "magenta") {col <- "\\e[1;35m"}
    else if (color == "cyan")    {col <- "\\e[1;36m"}
    else if (color == "white")   {col <- "\\e[1;37m"}
    end <- "\\e[0m"
    system2("echo", args = c("-e", paste0("\"", col, text, end, "\"")))
  } else {
    stop("Running on an unknown OS\n")
  }
}

# Compute threshold from cont distribution
threshold <- function(xSeries,
                      all_names,
                      adapt = threshold_adapt,
                      thr_val = threshold_value,
                      out_folder = out_dir)
{
  # Get series_ID and counter
  series_ID <- attr(xSeries, "own_name")
  item <- which(all_names == series_ID)
  position <- paste0("[", item, "/", length(all_names), "]")
  echo(paste("\nxSeries", position, series_ID), "yellow")
  
  # Adaptive expression threshold
  if (adapt == "true") {
    # Take the log2 of counts
    only_counts <- log2(countMatrix(xSeries)[,-1] + 1)
    # Make box-plots of count distributions
    r4tcpl::savePlots(
      \(){boxplot(only_counts)},
      figure_Name = paste0(series_ID, "_boxplot"),
      figure_Folder = file.path(out_folder, series_ID),
      pdf_out = FALSE)
    # Find the expression threshold adaptively
    # Filter the dataset by keeping only those genes that are detected in the
    # majority of the samples, compute their average expression, then use that
    # distribution of mean log-counts to fit the GMM.
    gmm <- r4tcpl::GMM_divide(
      rowMeans(only_counts)[rowSums(only_counts > 0) > N_selection(xSeries)/2],
      G = thr_val)
    # Set the new expression threshold as the right-most decision boundary
    gmm$boundary[thr_val*(thr_val-1)/2] |> unname() -> thr
    # Make density plots with GMM overlaid
    r4tcpl::savePlots(
      \(){
        # Density curves
        r4tcpl::count_density(only_counts,
                              remove_zeros = TRUE,
                              xlim = c(-1,10),
                              col = "gray20",
                              titles = c(paste0("Kernel Density Plot\n",
                                                series_ID), ""))
        # Plot the GMM
        for (i in 1:thr_val) {
          lines(gmm$x, gmm$components[,i], col="dodgerblue")
        }
        lines(gmm$x, rowSums(gmm$components), col="firebrick2")
        # Plot the expression threshold
        y_lim <- par("yaxp")[2]
        lines(c(thr, thr), c(0, 1.5*y_lim), col="darkslategray", lty="dashed")
        original_adj <- par("adj") # Store the original value of 'adj'
        par(adj = 0) # Set text justification to left
        text(x = thr + 0.3, y = 0.8*y_lim,
             labels = paste("Decision Boundary =", round(thr, digits = 2)),
             cex = 1.1)
        par(adj = original_adj) # Restore the original 'adj' value
      },
      figure_Name = paste0(series_ID, "_threshold"),
      figure_Folder = file.path(out_folder, series_ID))
    # Reject threshold if too low (to be conservative)
    if (thr < 1) {
      cat("\nWARNING:\n Adaptive threshold from GMM returned",
          round(thr, digits = 2), "...been coerced to 1.\n")
      thr <- 1
    } else {cat("\nAdaptive expression threshold set to: thr =", thr, "\n")}
  } else {
    # Fixed expression threshold (non-adaptive mode)
    thr <- thr_val
  }
  return(thr)
}




