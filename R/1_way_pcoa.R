#' Perform a PCoA analysis of one variable
#' @import ggplot2
#' @import stats
#' @importFrom vegan adonis
#' @importFrom ape pcoa
#' @importFrom qiimer dist_subset
#'
#' @param sample_mapping The sample table as a data frame.
#' @param distance_matrix A matrix of object type dist.
#'     The matrix will be trimmed to include only the samples in `s`.
#' @param var The primary variable of interest (string) to be plotted
#' @param strata If data is stratified, the name of the column containing the strata
#' @param ci Confidence interval to plot ellipse around. NULL if not plotted.
pcoa_1way <- function( sample_mapping, distance_matrix, var="StudyGroup", strata=NULL, ci=NULL ) {
  # Select columns of interest from distance matrix
  distances <- dist_subset( distance_matrix, sample_mapping$SampleID )

  res <- list()

  # Run the actual PCoA ordination
  res$pcoa <- pcoa( distances )

  # Save distances and sample info back into result
  res$distances <- distances
  res$df <- cbind(sample_mapping, res$pcoa$vectors[,1:2])
  test.pca.percExp <- round( res$pcoa$values[,1] / sum(res$pcoa$values[,1]) * 100, 2 )

  # Make a simple plot, save in result
  res$plot <- ggplot( res$df, aes_string( x="Axis.1", y="Axis.2", color=var )) +
    geom_point() +
    theme_classic() +
    xlab( sprintf("PC1 %s%s", test.pca.percExp[1], "%" )) +
    ylab(sprintf( "PC2 %s%s", test.pca.percExp[2], "%" ))
  # scale_colour_brewer(palette="Set1") +
  # geom_text(aes(label = SampleID )) +  ( add , label=s$SampleID to aes_string if using this)


  if (!is.null(ci)) {
    res$plot <- res$plot + stat_ellipse(level=ci, linetype=2)
  }

  # SPLIT FUNCTIONS INTO PLOT AND ADONIS !!!

  # Run an adonis test, save in result
  adonis_formula <- formula( paste("distances ~", var ))
  strata_arg <- if (is.null(strata)) NULL else res$df[[strata]]
  res$test <- adonis(adonis_formula, data=sample_mapping, strata=strata_arg)

  # Fix the test results:
  # Print the actual formula instead of "adonis_formula"
  # Replace the strata, too, for readability
  adonis_call <- as.list(res$test$call)
  adonis_call$formula <- adonis_formula
  # Have to do some R metaprogramming to make this look right...
  adonis_call$strata <- if (is.null(strata)) NULL else substitute(sample_mapping$strata, list(strata=strata))
  res$test$call <- as.call(adonis_call)

  class(res) <- "pcoa_1way"
  res
}

print.pcoa_1way <- function (x) {
  print(x$plot)
}
