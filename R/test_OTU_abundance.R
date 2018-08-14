#' Test OTUs for differences in abundance
#' @import ggplot2
#' @import stats
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange
#'
#' @param s A data frame of sample info.
#' @param cts A matrix of OTU counts (rows) for each sample (columns).
#' @param var A factor of sample groups to be tested.
#' @param a A chearacter vector of OTU assignments.
#' @param min_fraction Minimum fraction of samples that must contain counts of an OTU to be reported.
test_otu_abundance <- function (s, cts, var="study_group", a=NULL, min_fraction=0.5) {
  # Limit OTU table to samples in provided data frame
  cts <- cts[, s$SampleID]

  # Matrix of OTU proportions
  props <- sweep(cts, 2, colSums(cts), `/`)

  # Empty result
  res <- list()

  # Detect OTUs present in < 5 samples
  frac_present <- function (x) sum(x > 0) / length(x)
  too_rare <- apply(cts, 1, function (x) frac_present(x) < min_fraction)

  res$tests <- apply(cts[!too_rare,], 1, function (x) kruskal.test(x, s[[var]]))

  res$df <- data.frame(
    OtuID = names(res$tests),
    Taxon = a[names(res$tests)],
    Stat = sapply(res$tests, `[[`, "statistic"),
    Pval = sapply(res$tests, `[[`, "p.value"))

  # res$df <- within(res$df, {
  #   Fdr <- p.adjust(Pval, method="fdr")
  # })
  #this avoids appearing (unjustly) as a global variable to R CMD CHECK
  res$df$Fdr <- p.adjust(res$df$Pval, method="fdr")

  res$df <- arrange(res$df, Pval)

  plot_otus <- as.character(subset(res$df, Pval < 0.05)$OtuID)
  if (length(plot_otus) > 20) {
    plot_otus <- plot_otus[1:20]
  }

  plotdf <- melt(props[plot_otus,], varnames=c("OtuID", "SampleID"), value.name="Proportion")
  plotdf <- merge(plotdf, res$df, by="OtuID", all.x=T, all.y=F)
  plotdf <- merge(plotdf, s, by="SampleID", all.x=T, all.y=F)
# again avoiding an unjust NOTE
  # plotdf <- within(plotdf, {
  #   Label <- factor(paste(OtuID, Taxon))
  # })
  plotdf$Label <- factor(paste(OtuID, Taxon))

  plotdf$Label <- reorder(plotdf$Label, plotdf$Pval)

  res$plot <- ggplot(plotdf) +
    geom_boxplot(aes_string(x=var, y="Proportion")) +
    facet_wrap(~ Label, scales="free_y", ncol=2) +
    theme_classic() +
    theme(
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8)
    )

  class(res) <- "otu_abundance_tests"
  res
}
print.otu_abundance_tests <- function (x) {
  print(x$plot)
  df <- subset(x$df, Pval < 0.05)
  if (nrow(df) > 20) {
    df <- df[1:20,]
  }
  print(df)
}
