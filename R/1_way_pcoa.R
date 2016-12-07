require(proto)
require(vegan)
require(ape)
require(ggplot2)
require(plyr)
require(reshape2)


#' Perform a PCoA analysis of one variable
#' @param s The sample table
#' @param d A distance matrix.  The matrix will be trimmed to include only the samples in `s`.
#' @param var The variable of interest (string)
#' @param strata If data is stratified, the name of the column containing the strata
#' @param ci Confidence interval to plot ellipse around. NULL if not plotted.
pcoa_1way <- function(s, d, var="StudyGroup", strata=NULL, ci=NULL) {
  # Select columns of interest from distance matrix
  distances <- dist_subset(d, s$SampleID)
  
  res <- list()
  
  # Run the actual PCoA ordination
  res$pcoa <- pcoa(distances)
  
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


# Borrowed from Josef Fruehwald, "The Phonological Influence on Phonetic Change" - UPenn 2013
# StatEllipse <- proto(ggplot2:::Stat,
#   {
#     required_aes <- c("x", "y")
#     default_geom <- function(.) GeomPath
#     objname <- "ellipse"
#     
#     calculate_groups <- function(., data, scales, ...){
#       .super$calculate_groups(., data, scales,...)
#     }
#     calculate <- function(., data, scales, level = 0.75, segments = 51,...){
#       dfn <- 2
#       dfd <- length(data$x) - 1
#       if (dfd < 3){
#         ellipse <- rbind(c(NA,NA))	
#       } else {
#         require(MASS)
#         v <- cov.trob(cbind(data$x, data$y))
#         shape <- v$cov
#         center <- v$center
#         radius <- sqrt(dfn * qf(level, dfn, dfd))
#         angles <- (0:segments) * 2 * pi/segments
#         unit.circle <- cbind(cos(angles), sin(angles))
#         ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
#       }
#       
#       ellipse <- as.data.frame(ellipse)
#       colnames(ellipse) <- c("x","y")
#       return(ellipse)
#     }
#   }
# )
# 
# stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
#   StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
# }
