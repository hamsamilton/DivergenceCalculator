#' Standardize read depth across samples
#'
#' @description
#' Standardizes read depth across samples by downsampling to a specified number of reads.
#' Columns with fewer reads than the specified number are skipped.
#'
#' @param df A data frame or matrix where rows are genes and columns are samples. Row names should be gene identifiers.
#' @param numreads Number of reads to downsample to. Default is 1,000,000.
#' @param ntimes Number of times to perform the downsampling. Default is 1.
#'
#' @return A data frame with standardized read counts.
#'
#' @importFrom dplyr left_join select
#'
#' @examples
#' \dontrun{
#' # Standardize read depth to 500,000 reads
#' standardized_counts <- StandardizeReadDepth(expression_matrix, numreads = 5e5)
#'
#' # Take the average of 3 downsamples to 3e6 trsfd
#' standardized_counts <- StandardizeReadDepth(expression_matrix, numreads = 1e6, ntimes = 3)
#' }
#'
#' @export
StandardizeReadDepth <- function(df, numreads = 1e6, ntimes = 1, parallel = FALSE, cores = NULL, verbose = TRUE) {
  
  # Start timing the whole process
  total_start_time <- Sys.time()
  if (verbose) {
    cat(sprintf("Starting StandardizeReadDepth at %s\n", format(total_start_time, "%H:%M:%S")))
    cat(sprintf("Processing %d samples with %d iterations each\n", ncol(df), ntimes))
    cat(sprintf("Target read depth: %s reads per sample\n", format(numreads, big.mark=",")))
  }
    
    # Setup parallel processing
    if (is.null(cores)) {
      cores <- max(1, detectCores() - 1) # Use all cores minus 1 by default

    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    if (verbose) {
      cat(sprintf("Running in parallel mode with %d cores\n", cores))
    }
    
    on.exit({
      stopCluster(cl)
      if (verbose) {
        cat("Parallel cluster shut down successfully\n")
      }
    }) # Ensure cluster is stopped when function exits
  } else {

    if (verbose) {
      cat("Running in sequential mode\n")
    }
  }
  
  df <- as.data.frame(df, row.names = rownames(df))
  gene_symbols <- rownames(df)
  
  # Helper function that does one individual column
  standardizereaddepth <- function(rnacol, clname = "rsmpld", ntimes, numreads = 1e6, verbose = TRUE) {
    sample_start_time <- Sys.time()
    if (verbose) {
      cat(sprintf("[%s] Started processing sample: %s\n", format(sample_start_time, "%H:%M:%S"), clname))
    }
    
    hlddf <- data.frame(row = 1:length(rnacol))
    hlddf$row <- as.factor(hlddf$row)
    totreads <- sum(rnacol)
    
    # Check if there are enough reads to sample
    if (totreads < numreads) {
      warning(paste("Warning: Sample",
                    clname,
                    "has fewer reads (", totreads, ") than requested to sample (", numreads, "). Skipping this sample."))
      if (verbose) {
        cat(sprintf("[%s] Skipped sample %s (insufficient reads)\n", format(Sys.time(), "%H:%M:%S"), clname))
      }
      return(NULL)  # Return NULL instead of a dataframe with NAs
    }
    
    lklihood <- rnacol / totreads
    lklihood[is.na(lklihood)] <- 0
    
    if (verbose && ntimes > 1) {
      cat(sprintf("[%s] Sample %s: Preparing for %d iterations to be averaged\n", 
                  format(Sys.time(), "%H:%M:%S"), clname, ntimes))
    }
    
    # Initialize a matrix to store counts from each iteration
    # This will be used to calculate the average later
    iteration_counts <- matrix(0, nrow = length(rnacol), ncol = ntimes)
    successful_iterations <- 0
    
    for (i in 1:ntimes) {
      iter_start_time <- Sys.time()
      if (verbose && ntimes > 1) {
        cat(sprintf("[%s] Sample %s: Starting iteration %d/%d\n", 
                    format(iter_start_time, "%H:%M:%S"), clname, i, ntimes))
      }
      
      # Create a vector where each gene index is repeated by its count
      # This simulates the actual reads in the sample
      all_reads <- unlist(lapply(1:length(rnacol), function(idx) {
        # Ensure we don't create negative or NA values when rounding
        count <- max(0, round(rnacol[idx]))
        if (is.na(count) || count == 0) return(NULL)
        rep(idx, count)
      }))
      
      # Now sample from all actual reads without replacement
      if (length(all_reads) >= numreads) {
        if (verbose) {
          cat(sprintf("[%s] Sample %s: Sampling %s reads from %s available reads\n", 
                      format(Sys.time(), "%H:%M:%S"), clname, 
                      format(numreads, big.mark=","), 
                      format(length(all_reads), big.mark=",")))
        }
        
        sampled_reads <- sample(all_reads,
                                size = numreads,
                                replace = F)
        
        # Convert to table to count occurrences
        read_counts <- table(factor(sampled_reads, levels = 1:length(rnacol)))
        
        # Store the counts in our matrix
        iteration_counts[, i] <- as.numeric(read_counts)
        successful_iterations <- successful_iterations + 1
        
        iter_end_time <- Sys.time()
        if (verbose && ntimes > 1) {
          time_taken <- as.numeric(difftime(iter_end_time, iter_start_time, units="secs"))
          cat(sprintf("[%s] Sample %s: Completed iteration %d/%d (%.2f seconds)\n", 
                      format(iter_end_time, "%H:%M:%S"), clname, i, ntimes, time_taken))
        }
      } else {
        warning(paste("Warning: Sample", clname, "has fewer actual reads (", length(all_reads),
                      ") than requested to sample (", numreads, "). Skipping this iteration."))
        if (verbose) {
          cat(sprintf("[%s] Sample %s: Skipped iteration %d (insufficient reads)\n", 
                      format(Sys.time(), "%H:%M:%S"), clname, i))
        }
        # Skip this iteration without adding to the counts
        next
      }
    }
    
    # Calculate the average counts across all successful iterations
    if (successful_iterations > 0) {
      if (verbose && ntimes > 1) {
        cat(sprintf("[%s] Sample %s: Calculating average across %d successful iterations\n", 
                    format(Sys.time(), "%H:%M:%S"), clname, successful_iterations))
      }
      
      # Calculate average counts
      avg_counts <- rowSums(iteration_counts) / successful_iterations
      
      # Add average counts to the output dataframe
      avg_df <- data.frame(
        row = factor(1:length(rnacol)),
        count = avg_counts
      )
      colnames(avg_df)[2] <- clname
      
      # Join with the main dataframe
      hlddf <- dplyr::left_join(hlddf, avg_df, by = "row")
    } else {
      if (verbose) {
        cat(sprintf("[%s] Sample %s: No successful iterations completed\n", 
                    format(Sys.time(), "%H:%M:%S"), clname))
      }
      return(NULL)
    }
    
    # Replace NAs with 0s
    hlddf[is.na(hlddf)] <- 0
    
    sample_end_time <- Sys.time()
    time_taken <- as.numeric(difftime(sample_end_time, sample_start_time, units="secs"))
    if (verbose) {
      cat(sprintf("[%s] Completed sample %s in %.2f seconds\n", 
                  format(sample_end_time, "%H:%M:%S"), clname, time_taken))
    }
    
    return(hlddf)
  }
  
  # Initialize the dataframe where the downsampled reads are stored
  init_df <- data.frame(row = 1:nrow(df))
  init_df$row <- as.factor(init_df$row) # store indices to match reads back to genes
  
  # Calculate and display estimated time
  if (verbose) {
    cat("\n--- Processing samples ---\n")
  }
  
  # Timer for progress tracking
  start_time <- Sys.time()
  
  if (parallel) {
    # Process each column in parallel
    if (verbose) {
      cat(sprintf("Setting up parallel processing for %d samples\n", ncol(df)))
    }
    
    # Export the verbose parameter to the workers
    results <- foreach(i = 1:ncol(df), .packages = c("dplyr")) %dopar% {
      rnacol <- df[, i]
      clname <- colnames(df)[i]
      
      standardizereaddepth(
        rnacol = rnacol,
        clname = clname,
        ntimes = ntimes,
        numreads = numreads,
        verbose = verbose
      )
    }
    
    # Join all results with the initial dataframe
    if (verbose) {
      cat(sprintf("[%s] Combining results from all samples\n", format(Sys.time(), "%H:%M:%S")))
    }
    
    final_df <- init_df
    successful_samples <- 0
    for (i in 1:length(results)) {
      if (!is.null(results[[i]])) {
        final_df <- dplyr::left_join(final_df, results[[i]], by = "row")
        successful_samples <- successful_samples + 1
      }
    }
    
    if (verbose) {
      cat(sprintf("[%s] Successfully processed %d out of %d samples\n", 
                  format(Sys.time(), "%H:%M:%S"), successful_samples, ncol(df)))
    }
  } else {
    # Sequential processing of each column
    final_df <- init_df
    completed <- 0
    skipped <- 0
    
    for (i in 1:ncol(df)) {
      rnacol <- df[, i]
      clname <- colnames(df)[i]
      
      # Calculate and display progress
      if (verbose && i > 1) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
        avg_time_per_sample <- elapsed / (i - 1)
        remaining_samples <- ncol(df) - (i - 1)
        est_remaining_time <- avg_time_per_sample * remaining_samples
        
        cat(sprintf("\n[%s] Progress: %d/%d samples (%.1f%%) - Est. time remaining: %s\n", 
                    format(Sys.time(), "%H:%M:%S"),
                    i - 1, ncol(df), 
                    (i - 1) / ncol(df) * 100,
                    format(as.POSIXct(Sys.time()) + est_remaining_time, "%H:%M:%S")))
      }
      
      r50out <- standardizereaddepth(
        rnacol = rnacol,
        clname = clname,
        ntimes = ntimes,
        numreads = numreads,
        verbose = verbose
      )
      
      # Only join if r50out is not NULL (column had enough reads)
      if (!is.null(r50out)) {
        final_df <- dplyr::left_join(final_df, r50out, by = "row")
        completed <- completed + 1
      } else {
        skipped <- skipped + 1
      }
    }
    
    if (verbose) {
      cat(sprintf("\n[%s] Processing complete - %d samples processed, %d samples skipped\n", 
                  format(Sys.time(), "%H:%M:%S"), completed, skipped))
    }
  }
  
  # Set row names and clean up the dataframe
  rownames(final_df) <- gene_symbols
  final_df <- final_df %>% dplyr::select(-row)
  
  # Total execution time
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, total_start_time, units="mins"))
  
  if (verbose) {
    cat(sprintf("\n=== StandardizeReadDepth Complete ===\n"))
    cat(sprintf("Started: %s\n", format(total_start_time, "%H:%M:%S")))
    cat(sprintf("Finished: %s\n", format(end_time, "%H:%M:%S")))
    cat(sprintf("Total execution time: %.2f minutes\n", total_time))
    cat(sprintf("Final output dimensions: %d genes × %d samples\n", 
                nrow(final_df), ncol(final_df)))
  }
  
  return(final_df)
}
#' Create gene expression distribution plots
#'
#' @description
#' Creates histogram plots of gene expression distribution for each sample.
#'
#' @param GeneExpressionMat A data frame or matrix where rows are genes and columns are samples.
#'
#' @return A list of ggplot objects, one for each sample.
#'
#' @details
#' This function creates distribution plots for gene expression data. It removes
#' zero-count genes, transforms data to log2 scale, and adds reference lines for
#' median expression and defined thresholds (LowP50 and HighP50).
#'
#' Note: The variables LowP50 and HighP50 need to be defined before calling this function.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs coord_cartesian theme_classic theme element_text
#'
#' @examples
#' \dontrun{
#' # Define threshold values
#' LowP50 <- 2
#' HighP50 <- 8
#'
#' # Create distribution plots
#' plots <- mkGeneDistPlots(expression_matrix)
#'
#' # Save plots to a directory
#' dir.create("expression_plots", showWarnings = FALSE)
#' for (i in 1:length(plots)) {
#'   ggsave(paste0("expression_plots/", names(plots)[i], ".pdf"),
#'          plots[[i]], width = 8, height = 6)
#' }
#' }
#'
#' @export
mkGeneDistPlots = function (GeneExpressionMat) {
  SampleNames = colnames(GeneExpressionMat)
  
  
  # Define histogram breaks (same as used in the plots)
  bin_breaks = seq(0.1, 15, length.out = 41)
  
  # Calculate per-sample maximum histogram counts using only values in [0.1, 15]
  sample_max_counts = sapply(SampleNames, function(SampleName) {
    SampleVec = GeneExpressionMat[[SampleName]]
    SampleVec = SampleVec[SampleVec != 0]
    SampleVec = sort(SampleVec)
    SampleVec = log2(SampleVec + 1)
    # Filter to values that fall within the defined range
    SampleVec = SampleVec[SampleVec >= 0.1 & SampleVec <= 15]
    counts = hist(SampleVec, breaks = bin_breaks, plot = FALSE)$counts
    print(SampleName)
    max(counts)
  })
  global_y_max = max(sample_max_counts)
  
  P50Plots = lapply(X   = SampleNames,
                    FUN = function(SampleName) {

                      # Preprocess Gene Expression Vector
                      SampleVec = GeneExpressionMat[[SampleName]]
                      SampleVec = SampleVec[SampleVec != 0]
                      SampleVec = sort(SampleVec)

                      n = length(SampleVec) # How many genes were detected

                      # Transform into the log scale so that the distribution can be visualized
                      SampleVec = log2(SampleVec + 1)
                      MedianGeneExpression = median(SampleVec)
                      LowP50 <- mean(SampleVec[1:floor(n/2)])
                      HighP50 <- mean(SampleVec[(ceiling(n/2) + 1):n])

                      P50Plot = ggplot(data    = data.frame(Values = SampleVec),
                                       mapping = aes(x = Values)) +
                        geom_histogram(fill  = "skyblue",
                                       bins  = 40,
                                       color = "black") +
                        geom_vline(xintercept = MedianGeneExpression,
                                   color      = "black",
                                   linetype   = "dashed",
                                   linewidth  = 1) +
                        geom_vline(xintercept = HighP50,
                                   color      = "red",
                                   linetype   = "dashed",
                                   linewidth  = 1) +
                        geom_vline(xintercept = LowP50,
                                   color      = "blue",
                                   linetype   = "dashed",
                                   linewidth  = 1) +
                        labs(title = SampleName,
                             x     = "Gene Expression (log2)",
                             y     = "# of Genes") +
                        coord_cartesian(xlim   = c(0.1,15),
                                        ylim   = c(0, global_y_max + (global_y_max * .05)), 
                                        expand = F) +
                        theme_classic() +
                        theme(title      = element_text(size  = 20,
                                                        hjust = 0.5),
                              axis.title = element_text(size = 20),
                              axis.text  = element_text(size = 18))
                      # TODO Add a legend top right of plot.
                      return(P50Plot)
                    })
  return(P50Plots)
}

#' Calculate transcriptomic divergence metrics
#'
#' @description
#' Calculates various metrics of transcriptomic divergence across samples.
#'
#' @param cpmdf A depth-corrected data frame (e.g., CPM or FPKM), where rows are genes and columns are samples.
#' @param metric The divergence metric to calculate. Default is "P50/P50". Options include:
#'   \itemize{
#'     \item "P50/P50": Ratio of sum of top 50% expressed genes to bottom 50%
#'     \item "P20": Ratio of sum of top 20% expressed genes to bottom 20%
#'     \item "P10": Ratio of sum of top 10% expressed genes to bottom 10%
#'     \item "sd": Standard deviation of gene expression values
#'     \item "kurtosis": Kurtosis of gene expression distribution
#'     \item "gini": Gini coefficient of gene expression values
#'   }
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item dvrgnc: The calculated divergence value for each sample
#'     \item clnames: The sample names
#'   }
#'
#' @details
#' This function automatically removes genes with zero counts before calculating metrics.
#'
#' @importFrom e1071 kurtosis
#' @importFrom reldist gini
#'
#' @examples
#' \dontrun{
#' # Calculate P50/P50 (default)
#' p50_divergence <- mkdvrgnc(normalized_counts)
#'
#' # Calculate Gini coefficient
#' gini_divergence <- mkdvrgnc(normalized_counts, metric = "gini")
#'
#' # Plot divergence metrics
#' library(ggplot2)
#' ggplot(p50_divergence, aes(x = clnames, y = dvrgnc)) +
#'   geom_bar(stat = "identity") +
#'   theme_minimal() +
#'   labs(title = "P50/P50 Divergence by Sample",
#'        x = "Sample",
#'        y = "P50/P50 Ratio")
#' }
#'
#' @export
CalcDivergence <- function(cpmdf, metric = "P50/P50") {

  dvrgnc <- c() # Init where divergence values are stored
  clnames <- c() # init where column names are stored
  for (i in 1:ncol(cpmdf)){
    col = cpmdf[,i]
    clnames <- append(clnames,colnames(cpmdf)[i])
    col = col[which(col > 0)] # get rid of 0 count

    #calcp50/p50
    if(metric == "P50/P50"){
      col = col[order(col)]
      midindex = floor(length(col)*.5)
      coltop = col[1:midindex]
      colbot = col[(midindex+1):length(col)]
      dvrgnc = append(dvrgnc,sum(colbot) / sum(coltop))}

    if(metric == "P20"){
      col = col[order(col)]
      botindex = floor(length(col)*.2)
      topindex = length(col) - botindex
      coltop = col[1:botindex]
      colbot = col[topindex:length(col)]
      dvrgnc = append(dvrgnc,sum(colbot) / sum(coltop))}

    if(metric == "P10"){
      col = col[order(col)]
      botindex = floor(length(col)*.1)
      topindex = length(col) - botindex
      coltop = col[1:botindex]
      colbot = col[topindex:length(col)]
      dvrgnc = append(dvrgnc,sum(colbot) / sum(coltop))}


    if(metric == "sd")       dvrgnc <- append(dvrgnc,sd(col))
    if(metric == "kurtosis") dvrgnc <- append(dvrgnc,e1071::kurtosis(col))
    if(metric == "gini")     dvrgnc <- append(dvrgnc,reldist::gini(col))
  }
  dvrgncframe <- data.frame(dvrgnc,clnames)
  return(dvrgncframe)
}
