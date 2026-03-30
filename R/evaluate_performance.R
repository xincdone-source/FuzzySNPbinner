#' Evaluate Algorithmic Performance (Simple Summary)
#'
#' 用于快速比较：模糊版 vs 原版 的 bin 数量、平均BIN大小等。
#' @param bins_csv 模糊版输出的 bins.csv（或目录下汇总）
#' @param baseline_bins_csv 原版SNPbinner的 bins.csv（可选）
#' @export
evaluate_performance <- function(bins_csv, baseline_bins_csv = NULL) {
  if (!file.exists(bins_csv)) stop("File not found: ", bins_csv)
  df <- readr::read_csv(bins_csv, show_col_types = FALSE)
  msg <- list(
    fuzzy_bin_count = nrow(df),
    fuzzy_bin_size_mean = mean(df$Bin_Size, na.rm = TRUE),
    fuzzy_bin_size_median = median(df$Bin_Size, na.rm = TRUE)
  )
  if (!is.null(baseline_bins_csv) && file.exists(baseline_bins_csv)) {
    df0 <- readr::read_csv(baseline_bins_csv, show_col_types = FALSE)
    msg$baseline_bin_count <- nrow(df0)
    msg$baseline_bin_size_mean <- mean(df0$Bin_Size, na.rm = TRUE)
    msg$baseline_bin_size_median <- median(df0$Bin_Size, na.rm = TRUE)
  }
  return(msg)
}
