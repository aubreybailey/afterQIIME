#' Perform a PCoA analysis of one variable
#' @param cts The observation counts matrix
#' @param mapping_subset The data frame or subset thereof containing the samples to be analyzed.
#' @param var The variable of interest (string) (can be a vector)
#' @param strata If data is stratified, the name of the column containing the strata
#' @param ci Confidence interval to plot ellipse around. NULL if not plotted.
rarefaction_diversity <- function(cts=cts, mapping_subset, var="StudyGroup", rarefaction_depth=1000, do_test=TRUE) {
  res <- list()
  rare_for_test <- rarefy(t(cts), rarefaction_depth)
  mapping_subset$Rare <- rare_for_test[mapping_subset$SampleID]
  
  res$plot <- ggplot(mapping_subset) +
    geom_boxplot(aes_string(x=var,y="Rare") ) +
    labs(y="Observed OTUs",
         x=var,
         title=paste("Number of unique OTUs rarefied to", rarefaction_depth, "observations per sample" )) +
    theme_classic()+
    theme(axis.text.x = element_text(angle=60, hjust=1))
   
  if(do_test){
      res$test <- t.test(as.formula(paste("Rare ~", paste(var, collapse="+"))), data=mapping_subset)
  }

  class(res) <- "rarefaction_diversity"
  res
}


