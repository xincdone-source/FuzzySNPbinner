#' Preview Fuzzy HMM Parameters for a Given TSV
#'
#' 预览基于模糊数学自动推断的隐马尔可夫模型(HMM)参数
#' @param tsv 输入TSV
#' @param preview_k 窗口bp，默认200000
#' @return 包含 HMM 基础参数和覆盖参数的列表
#' @export
plot_fuzzy_params <- function(tsv, preview_k = 200000L) {
  if (!file.exists(tsv)) stop("File not found: ", tsv)
  fm <- system.file("python/fuzzy_module.py", package = "FuzzySNPbinner")
  su <- system.file("python/stats_utils.py", package = "FuzzySNPbinner")
  reticulate::source_python(fm)
  reticulate::source_python(su)

  stats <- estimate_stats(tsv, window_bp = as.integer(preview_k))
  mu <- calibrate_memberships(stats)

  base_hmm <- infer_hmm_base(mu, stats)
  overrides <- infer_hmm_overrides(mu, stats)

  message("=== 模糊 HMM 参数自动推断结果 ===")
  message(sprintf("HMM 基准参数: 预期同型率(homo)=%.3f, 预期交叉点(cross_count)=%d",
                  base_hmm$predicted_homogeneity, base_hmm$predicted_cross_count))
  message(sprintf("HMM 修正参数: 发射概率(intra_p)=%.3f, 转移惩罚乘数(trans_mult)=%.3f",
                  overrides$intra_p, overrides$trans_mult))

  invisible(list(base = base_hmm, overrides = overrides))
}
