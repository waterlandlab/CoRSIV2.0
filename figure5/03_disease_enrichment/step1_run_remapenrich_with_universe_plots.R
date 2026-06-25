suppressPackageStartupMessages({
  library(ReMapEnrich)
  library(data.table)
  library(GenomicRanges)
})

query_files <- c(
  "./query/Black_CpGge5_Rangegt20.bed",
  "./query/White_CpGge5_Rangegt20.bed"
)

universe_files <- c(
  "./background/sorted.matched.controls.bed",
  "./background/random_background.bed"
)

catalog_bed <- "./remap2022_all_macs2_hg38_v1_0.bed"
out_dir <- "remapenrich_results_2022_with_universe"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(catalog_bed) || file.info(catalog_bed)$size == 0) {
  stop("Catalog file not found or empty: ", catalog_bed)
}

cat("Using catalog:", catalog_bed, "\n")
catalog_gr <- bedToGranges(catalog_bed)
cat("Catalog regions:", length(catalog_gr), "\n")

summary_rows <- list()

filter_query_for_universe <- function(query_gr, universe_gr) {
  q_width <- width(query_gr)
  u_width <- width(universe_gr)
  q_chr <- as.character(seqnames(query_gr))
  u_chr <- as.character(seqnames(universe_gr))

  shared_chr <- base::intersect(unique(q_chr), unique(u_chr))
  keep_u_chr <- u_chr %in% shared_chr
  keep_q_chr <- q_chr %in% shared_chr

  keep_u <- keep_u_chr & !is.na(u_width) & u_width > 0
  universe_gr <- universe_gr[keep_u]
  u_width <- u_width[keep_u]

  if (length(universe_gr) == 0) {
    return(list(
      query_filtered = query_gr[0],
      universe_filtered = universe_gr,
      n_query_removed = length(query_gr),
      n_universe_removed = sum(!keep_u),
      n_query_removed_chr = sum(!keep_q_chr),
      n_universe_removed_chr = sum(!keep_u_chr)
    ))
  }

  max_u <- max(u_width, na.rm = TRUE)

  keep_q <- keep_q_chr & !is.na(q_width) & q_width > 0 & q_width <= max_u
  query_filtered <- query_gr[keep_q]

  return(list(
    query_filtered = query_filtered,
    universe_filtered = universe_gr,
    n_query_removed = sum(!keep_q),
    n_universe_removed = sum(!keep_u),
    n_query_removed_chr = sum(!keep_q_chr),
    n_universe_removed_chr = sum(!keep_u_chr)
  ))
}

for (qf in query_files) {
  if (!file.exists(qf) || file.info(qf)$size == 0) {
    stop("Query file not found or empty: ", qf)
  }

  query_gr <- bedToGranges(qf)
  qbase <- sub("\\.bed$", "", basename(qf), ignore.case = TRUE)
  cat("\nQuery:", qf, "regions=", length(query_gr), "\n")

  for (uf in universe_files) {
    if (!file.exists(uf) || file.info(uf)$size == 0) {
      stop("Universe file not found or empty: ", uf)
    }

    universe_gr <- bedToGranges(uf)
    ubase <- sub("\\.bed$", "", basename(uf), ignore.case = TRUE)

    prep <- filter_query_for_universe(query_gr, universe_gr)
    query_use <- prep$query_filtered
    universe_use <- prep$universe_filtered

    run_id <- paste(qbase, "vs", ubase, sep = "__")
    out_csv <- file.path(out_dir, paste0(run_id, "_enrichment.csv"))
    top_csv <- file.path(out_dir, paste0(run_id, "_top20_by_qvalue.csv"))
    dotplot_pdf <- file.path(out_dir, paste0(run_id, "_dotplot.pdf"))

    cat("  Universe:", uf, "regions=", length(universe_gr), "\n")
    cat("    Width filter removed query=", prep$n_query_removed,
        " universe=", prep$n_universe_removed,
        " -> using query=", length(query_use),
        " universe=", length(universe_use), "\n")
    cat("    Chromosome filter removed query=", prep$n_query_removed_chr,
        " universe=", prep$n_universe_removed_chr, "\n")

    if (length(query_use) == 0 || length(universe_use) == 0) {
      warning("Skipping run due to zero regions after filtering: ", run_id)
      next
    }

    if (file.exists(out_csv) && file.info(out_csv)$size > 0) {
      cat("    Reusing existing enrichment:", out_csv, "\n")
      enrich_dt <- fread(out_csv)
    } else {
      enrich_df <- enrichment(
        query = query_use,
        catalog = catalog_gr,
        universe = universe_use,
        byChrom = FALSE,
        nCores = 1
      )
      enrich_dt <- as.data.table(enrich_df)
      fwrite(enrich_dt, out_csv)
    }
    if (!("q.value" %in% colnames(enrich_dt))) {
      stop("Expected column 'q.value' not found for run: ", run_id)
    }

    top20 <- enrich_dt[order(get("q.value"))][1:min(20, .N), ]
    fwrite(top20, top_csv)

    grDevices::pdf(dotplot_pdf, width = 10, height = 8)
    tryCatch({
      enrichmentDotPlot(as.data.frame(enrich_dt), top = 20, main = run_id)
    }, finally = {
      dev.off()
    })

    sig <- enrich_dt[!is.na(get("q.value")) & get("q.value") < 0.05, ]
    n_sig <- nrow(sig)
    min_q <- if (n_sig > 0) min(sig[["q.value"]], na.rm = TRUE) else NA_real_

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      query_file = qf,
      universe_file = uf,
      n_query_regions = length(query_gr),
      n_query_regions_used = length(query_use),
      n_query_regions_removed_by_width = prep$n_query_removed,
      n_query_regions_removed_by_chr = prep$n_query_removed_chr,
      n_universe_regions = length(universe_gr),
      n_universe_regions_used = length(universe_use),
      n_universe_regions_removed_invalid = prep$n_universe_removed,
      n_universe_regions_removed_by_chr = prep$n_universe_removed_chr,
      n_catalog_regions = length(catalog_gr),
      n_enrichment_rows = nrow(enrich_dt),
      n_significant_qvalue_lt_0_05 = n_sig,
      min_significant_qvalue = min_q,
      enrichment_csv = out_csv,
      top20_csv = top_csv,
      dotplot_pdf = dotplot_pdf,
      stringsAsFactors = FALSE
    )

    cat("    Wrote:", out_csv, "\n")
    cat("    Wrote:", dotplot_pdf, "\n")
  }
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
summary_csv <- file.path(out_dir, "run_summary.csv")
fwrite(summary_df, summary_csv)
cat("\nRun summary:", summary_csv, "\n")
print(summary_df)
