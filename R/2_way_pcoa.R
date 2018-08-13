#require(proto)
require(vegan)
require(ape)
require(ggplot2)
#require(plyr)
require(reshape2)

#' Perform a PCoA analysis of one variable
#' @param s The sample table
#' @param d A distance matrix.  The matrix will be trimmed to include only the samples in `s`.
#' @param var The variable of interest (string)
#' @param strata If data is stratified, the name of the column containing the strata
#' @param ci Confidence interval to plot ellipse around. NULL if not plotted.
#' uu_op is a dist object
#' s_op is a data frame
pcoa_1way <- function(s, d, color="study_group", shape="study_day", strata=NULL, ci=NULL) {
  uu_op <- d
  s_op <- s

  #They are guaranteed to have samples in the same order if you do this:
  uu_op <- dist_subset(uu_op, s_op$SampleID)

  pc <- pcoa(uu_op)
  pctvar <- round(pc$values$Relative_eig * 100)
  pcdf <- cbind(s_op, pc$vectors[s_op$SampleID,1:3])

  ggplot(pcdf) +
    geom_point(aes(x=Axis.1, y=Axis.2, color=DayOfTreatment, shape=study_group),  size=4) +
    scale_color_brewer(palette="Set1") +
    scale_shape_discrete(solid=F) +
    labs(
      x=paste0("PCoA axis 1 (", pctvar[1], "%)"),
      y=paste0("PCoA axis 2 (", pctvar[2], "%)"),
      color=s_op$color, shape=s_op$shape) +
    theme_classic()
    ##################################

  res <- list()

  # Save distances and sample info back into result
  res$distances <- distances
  res$df <- cbind(s, res$pcoa$vectors[,1:2])
  test.pca.percExp<-round(res$pcoa$values[,1]/sum(res$pcoa$values[,1])*100, 2)

  # Make a simple plot, save in result
  res$plot <- ggplot(res$df, aes_string(x="Axis.1", y="Axis.2", color=var)) +
    geom_point() + theme_classic() + xlab(sprintf("PC1 %s%s", test.pca.percExp[1], "%")) +
    ylab(sprintf("PC2 %s%s", test.pca.percExp[2], "%"))
  # scale_colour_brewer(palette="Set1") +
  # geom_text(aes(label = SampleID )) +  ( add , label=s$SampleID to aes_string if using this)


  if (!is.null(ci)) {
    res$plot <- res$plot + stat_ellipse(level=ci, linetype=2)
  }

  # Run an adonis test, save in result
  adonis_formula <- formula(paste("distances ~", var))
  strata_arg <- if (is.null(strata)) NULL else res$df[[strata]]
  res$test <- adonis(adonis_formula, data=s, strata=strata_arg)

  # Fix the test results:
  # Print the actual formula instead of "adonis_formula"
  # Replace the strata, too, for readability
  adonis_call <- as.list(res$test$call)
  adonis_call$formula <- adonis_formula
  # Have to do some R metaprogramming to make this look right...
  adonis_call$strata <- if (is.null(strata)) NULL else substitute(s$strata, list(strata=strata))
  res$test$call <- as.call(adonis_call)

  class(res) <- "pcoa_1way"
  res
}

print.pcoa_1way <- function (x) {
  print(x$plot)
}
