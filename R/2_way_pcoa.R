#' Perform a PCoA analysis of two variables
#' @import ggplot2
#' @import stats
#' @importFrom vegan adonis
#' @importFrom ape pcoa
#' @importFrom qiimer dist_subset
#'
#' @param sample_mapping The sample table as a data frame.
#' @param distance_matrix A matrix of object type dist.
#'     The matrix will be trimmed to include only the samples in `s`.
#' @param color The primary variable to be plotted
#' @param shape The secondary variable to be plotted.
#' @param strata If data is stratified, the name of the column containing the strata
#' @param ci Confidence interval to plot ellipse around. NULL if not plotted.
pcoa_2way <- function(sample_mapping, distance_matrix, color="study_group", shape="study_day", strata=NULL, ci=NULL) {

  #They are guaranteed to have samples in the same order if you do this:
  distances <- dist_subset(distance_matrix, sample_mapping$SampleID)

  pc <- pcoa(distance_matrix)
  pctvar <- round(pc$values$Relative_eig * 100)
  pcdf <- cbind(sample_mapping, pc$vectors[sample_mapping$SampleID,1:3])

  ggplot(pcdf) +
    geom_point(aes(x=Axis.1, y=Axis.2, color=DayOfTreatment, shape=study_group),  size=4) +
    scale_color_brewer(palette="Set1") +
    scale_shape_discrete(solid=F) +
    labs(
      x=paste0("PCoA axis 1 (", pctvar[1], "%)"),
      y=paste0("PCoA axis 2 (", pctvar[2], "%)"),
      color=sample_mapping$color, shape=sample_mapping$shape) +
    theme_classic()
    ##################################

  res <- list()

  # Save distances and sample info back into result
  res$distances <- distances
  res$df <- cbind(sample_mapping, res$pcoa$vectors[,1:2])
  test.pca.percExp<-round(res$pcoa$values[,1]/sum(res$pcoa$values[,1])*100, 2)

  # Make a simple plot, save in result
  res$plot <- ggplot(res$df, aes_string(x="Axis.1", y="Axis.2", color=var)) +
    geom_point() + theme_classic() + xlab(sprintf("PC1 %s%s", test.pca.percExp[1], "%")) +
    ylab(sprintf("PC2 %s%s", test.pca.percExp[2], "%"))
  # scale_colour_brewer(palette="Set1") +
  # geom_text(aes(label = SampleID )) +  ( add , label=sample_mapping$SampleID to aes_string if using this)


  if (!is.null(ci)) {
    res$plot <- res$plot + stat_ellipse(level=ci, linetype=2)
  }

  # Run an adonis test, save in result
  adonis_formula <- formula(paste("distances ~", var))
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
