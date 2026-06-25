suppressPackageStartupMessages({
  library(data.table)
})

input_pattern <- "remapenrich_results_2022_with_universe/*_enrichment.csv"
files <- Sys.glob(input_pattern)
if (length(files) != 4) {
  stop("Expected 4 enrichment files, found ", length(files), ". Pattern: ", input_pattern)
}

parse_run_label <- function(path) {
  nm <- basename(path)
  nm <- sub("_enrichment\\.csv$", "", nm)
  nm <- gsub("__vs__", " vs ", nm)
  nm <- gsub("_", " ", nm)
  nm <- gsub("CpGge5 Rangegt20", "CpG>=5 Range>20", nm)
  nm
}

lst <- lapply(files, function(f) {
  dt <- fread(f, select = c("category", "q.value"))
  dt[, run_label := parse_run_label(f)]
  dt
})

all_dt <- rbindlist(lst, fill = TRUE)
all_dt <- all_dt[is.finite(q.value) & q.value > 0]
all_dt[, neglog10_q := -log10(q.value)]

# One value per category/run used by the heatmap
heat_long <- all_dt[, .(neglog10_q = max(neglog10_q, na.rm = TRUE), q.value = min(q.value, na.rm = TRUE)), by = .(category, run_label)]

heat_wide <- dcast(heat_long, category ~ run_label, value.var = "neglog10_q")
run_cols <- setdiff(colnames(heat_wide), "category")
heat_wide[, max_neglog10_q := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = run_cols]
setorder(heat_wide, -max_neglog10_q)

out_long <- "remapenrich_results_2022_with_universe/all_results_overlap_heatmap_data_long.csv"
out_wide <- "remapenrich_results_2022_with_universe/all_results_overlap_heatmap_data_wide.csv"

fwrite(heat_long, out_long)
fwrite(heat_wide, out_wide)

cat("Saved:\n", out_long, "\n", out_wide, "\n", sep = "")
cat("Rows (unique categories): ", nrow(heat_wide), "\n", sep = "")
