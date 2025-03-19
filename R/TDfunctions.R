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
StandardizeReadDepth <- function(df, numreads = 1e6, ntimes = 1) {
  df <- as.data.frame(df,
                      row.names = rownames(df))
  gene_symbols = rownames(df)
  # initialize the dataframe where the downsampled reads are stored
  hlddf <- data.frame(row = 1:nrow(df))
  hlddf$row <- as.factor(hlddf$row) # store indices to match reads back to genes

  #helper function that does one individual column
  standardizereaddepth <- function(rnacol, clname = "rsmpld", ntimes, numreads = 1e6){

    hlddf <- data.frame(row = 1:length(rnacol))
    hlddf$row <- as.factor(hlddf$row)
    totreads <- sum(rnacol)

    # Check if there are enough reads to sample
    if(totreads < numreads){
      warning(paste("Warning: Sample",
                    clname,
                    "has fewer reads (", totreads, ") than requested to sample (", numreads, "). Skipping this sample."))
      return(NULL)  # Return NULL instead of a dataframe with NAs
    }

    lklihood <- rnacol / totreads
    lklihood[is.na(lklihood)] <- 0
    for (i in 1:ntimes){
      print(clname)

      # Create a vector where each gene index is repeated by its count
      # This simulates the actual reads in the sample
      all_reads <- unlist(lapply(1:length(rnacol), function(idx) {
        rep(idx, round(rnacol[idx]))
      }))

      # Now sample from all actual reads without replacement
      if(length(all_reads) >= numreads) {
        sampled_reads <- sample(all_reads,
                                size = numreads,
                                replace = FALSE)

        # Convert to table to count occurrences
        smpld <- table(sampled_reads) %>%
          as.data.frame()
        colnames(smpld) <- c("row", clname)

        # Ensure row is a factor to match hlddf
        smpld$row <- as.factor(smpld$row)

        hlddf <- left_join(hlddf, smpld, by = c("row"))
      } else {
        warning(paste("Warning: Sample", clname, "has fewer actual reads (", length(all_reads),
                      ") than requested to sample (", numreads, "). Skipping this iteration."))
        # Skip this iteration without adding the column
        next
      }
    }
    hlddf[is.na(hlddf)] <- 0
    return(hlddf)
  }

  # Process each column
  for(i in 1:ncol(df)){
    rnacol <- df[,i]
    clname = colnames(df)[i]
    r50out <- standardizereaddepth(rnacol   = rnacol,
                                   clname   = clname,
                                   ntimes   = ntimes,
                                   numreads = numreads)

    # Only join if r50out is not NULL (column had enough reads)
    if(!is.null(r50out)) {
      hlddf <- hlddf %>%
        left_join(r50out, by = "row")
    }
  }

  rownames(hlddf) <- gene_symbols
  hlddf <- hlddf %>% dplyr::select(-row)

  return(hlddf)
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
                      LowP50 = log2(LowP50 + 1)
                      HighP50 = log2(HighP50 + 1)

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
                                        ylim   = c(0, 2000), # TODO Auto standardize value across plots
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
