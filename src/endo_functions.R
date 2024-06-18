

# Colored echo for on-screen log in Bash CLI
echo <- function(text, color = "white") {
  # Check the operating system
  os_type <- Sys.info()["sysname"]
  if (os_type == "Windows") {
    warning("Colored echo is just for Bash!")
    cat(text, "\n")
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

# Borrowed from SeqLoader
# Any named list knows the names of all the elements it contains (under its
# 'names' attribute), but it doesn't know its own (even when it has one, e.g.,
# because it is itself an element of a named list)! So, this function sets as
# attribute for each element of a named list its own name (to access it later).
set_own_names <- function(parent_list) {
  names(parent_list) |> lapply(function(element_name) {
    attr(parent_list[[element_name]], "own_name") <- element_name
    return(parent_list[[element_name]])
  }) |> setNames(names(parent_list)) # Also keep original names in parent list
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

# 'gois_stats' is supposed to be a data.frame with:
#   - descriptive statistics about gene expression in a Series;
#   - gene annotation columns featuring (at least) the 'SYMBOL' key;
#   - the "own_name" character attribute of the associated Series;
#   - the "metadata" data frame attribute about Runs participating in the stats.
plot_barChart <- function(gois_stats,
                          family_name,
                          y_limit,
                          border = FALSE,
                          thr,
                          out_folder = out_dir)
{
  # Colors
  line_color <- "gray17"
  err_color <- "gray17"
  # Get series_ID, sample size, and average megareads per Run
  gois_stats |> attr("own_name") -> series_ID
  gois_stats |> attr("metadata") -> meta_table
  meta_table |> nrow() -> n
  if (meta_table |> hasName("read_count")) {
    meta_table |> dplyr::select(read_count) |>
      (\(z){colMeans(z)/1e6})() |> round(digits = 2) -> megareads
  } else {megareads <- NA}
  if (meta_table |> hasName("library_layout")) {
    meta_table$library_layout[1] -> libLay
  } else {libLay <- NA}
  # Log Messages
  cat("xSeries ", series_ID, "... ", sep = "")
  
  # Prepare the Frame
  gg_frame <-
    ggplot(data = gois_stats,
           aes(x = SYMBOL, y = Mean, fill = SYMBOL)) +
    theme_bw(base_size = 15, base_rect_size = 1.5) +
    theme(axis.text.x = element_text(size = 10, angle = 90,
                                     vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.position = "none") +
    scale_y_continuous(limits=c(0, y_limit),
                       expand = c(0.01, 0.02),
                       breaks = seq(0, y_limit, 0.5)) +
    xlab("Genes of Interest") +
    ylab(substitute(log[2]*(x+1), list(x = "TPM"))) +
    ggtitle(label = paste0(series_ID, " (n = ", n,
                           ", depth = ", megareads, " MegaReads, ",
                           libLay, " library) - GOI subset: ", family_name))
  # Draw the Bars
  if (border) {
    gg_bars <- geom_bar(stat = "identity", width = 0.7,
                        color = line_color, linewidth = 0.1)
  } else {
    gg_bars <- geom_bar(stat = "identity", width = 0.75)
  }
  gg_errorbar <- gg_frame + gg_bars +
    geom_errorbar(aes(ymin = Mean - Std_Dev,
                      ymax = Mean + Std_Dev),
                  linewidth = 1.0, width = 0.5, color = err_color)
  # Add the Expression Threshold
  gg_thr <- gg_errorbar +
    geom_hline(yintercept = thr[[series_ID]],
               linetype = "dashed",
               color = line_color,
               linewidth = 1)
  # Save the Chart
  r4tcpl::savePlots(
    \(){print(gg_thr)},
    width_px = 2000,
    figure_Name = paste0(series_ID, "_", family_name, "_chart"),
    figure_Folder = file.path(out_folder, series_ID))
  
  cat("done\n")
}

# Regenerate annotation from scratch
add_annotation <- function(gene_matrix, OrgDb_key = "ENSEMBL") {
  # See columns(org.Hs.eg.db) or keytypes(org.Hs.eg.db) for a complete list of
  # all possible annotations.
  org_db <- org.Hs.eg.db::org.Hs.eg.db
  annots <- AnnotationDbi::select(org_db,
                                  keys = gene_matrix[,"IDs"],
                                  columns = c("SYMBOL", "GENENAME", "GENETYPE"),
                                  keytype = OrgDb_key)
  # Warning: 'select()' returned 1:many mapping between keys and columns
  # ========>
  # Collapse the duplicated entries in the ID column and concatenate the
  # (unique) values in the remaining columns using a comma as a separator to
  # prevent rows from being added in the following join step.
  if (anyDuplicated(annots[,OrgDb_key])) {
    cat("\nWARNING:\n Multiple annotation entries corresponding to a single\n",
        OrgDb_key, "ID will be collapsed by a comma separator.\n")
    annots <- aggregate(. ~ get(OrgDb_key),
                        data = annots,
                        FUN = \(x)paste(unique(x), collapse = ","),
                        na.action = NULL)[,-1]
  }
  colnames(annots)[colnames(annots) == OrgDb_key] <- "IDs" 
  gene_matrix <- merge(annots, gene_matrix, by = "IDs", all.y = TRUE)
}








