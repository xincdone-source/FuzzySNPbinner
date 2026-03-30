#' Run fuzzy binning on crosspoints output
#' @param input_path  Path to .crosp.csv (output of crosspoints)
#' @param output_path Path to save .bin.csv
#' @param min_bin_size Minimum bin size in bp (default 500,000)
#' @export
fuzzy_run_bins <- function(input_path, output_path, min_bin_size = 500000) {
  py <- system.file("python", "fuzzy_bins.py", package = "FuzzySNPbinner")
  if (py == "") stop("fuzzy_bins.py not found in package")
  reticulate::source_python(py)

  # 确保是整数传入 Python
  bins_fuzzy(
    input_path  = input_path,
    output_path = output_path,
    min_bin_size = as.integer(min_bin_size)
  )
}
