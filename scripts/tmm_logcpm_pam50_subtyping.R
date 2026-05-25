if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("edgeR", "genefu")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

invisible(capture.output(
  suppressPackageStartupMessages(suppressWarnings({
    library(edgeR)
    library(genefu)
  }))
))


tmm_logcpm_tsv <- function(input_tsv, output_tsv, prior_count = 1.0, return_df = TRUE, save = TRUE) {
  # input: TSV with rownames in first column, donors x genes (counts)
  # output: TSV with rownames in first column, donors x genes (TMM + logCPM)
  
  cat("[tmm_logcpm_tsv] reading:", input_tsv, "\n"); flush.console()
  t0 <- proc.time()

  counts_df <- read.table(
    input_tsv,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    row.names = 1
  )
  
  cat("[tmm_logcpm_tsv] running edgeR DGEList + TMM + cpm(log=TRUE)\n"); flush.console()

  # edgeR expects genes x samples
  d <- DGEList(counts = t(as.matrix(counts_df)))
  d <- normLibSizes(d, method = "TMM")

  lcpm <- cpm(d, log = TRUE, prior.count = prior_count)

  # back to donors x genes
  out_df <- as.data.frame(t(lcpm))

  if (save == TRUE) {
    cat("[tmm_logcpm_tsv] writing:", output_tsv, "\n"); flush.console()
    write.table(out_df, file = output_tsv, sep = "\t", quote = FALSE, col.names = NA)
  }
  
  cat("[tmm_logcpm_tsv] done in sec:", (proc.time() - t0)[["elapsed"]], "\n"); flush.console()

  if (return_df == TRUE){
    return(out_df)
  }
}


pam50_genefu_subtyping <- function(expr_df) {

  expr <- as.matrix(expr_df)  # samples x genes

  data("pam50.robust", package="genefu")
  res <- molecular.subtyping(
    sbt.model  = "pam50",
    data       = expr,
    annot      = NULL,
    do.mapping = FALSE,
    verbose    = TRUE
  )

  write.table(res, file = "../data/genefu_result.tsv", sep = "\t", quote = FALSE, col.names = NA)

  subtype_df <- data.frame(subtype = as.character(res$subtype),
                           row.names = names(res$subtype))
  subtype_df$subtype[subtype_df$subtype == "Her2"] <- "ERBB2"

  return(subtype_df)
}
