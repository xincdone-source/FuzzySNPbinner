#' Run fuzzy crosspoints on crosspoints data with manual min_length and min_bin_size
#' @param input_path Path to the input .tsv file
#' @param output_path Path to the output .crosp.csv file
#' @param min_length Manually set minimum length (default = 5000)
#' @param min_bin_size Manually set minimum bin size (default = 500000)
#' @param seq_depth Average sequencing depth (e.g., 10 for 10X). Affects HMM confidence. Default 10.
#' @export
fuzzy_run_crosspoints <- function(input_path, output_path, min_length = 5000, min_bin_size = 500000, seq_depth = 10) {
  py <- system.file("python", "fuzzy_crosspoints.py", package = "FuzzySNPbinner")
  if (py == "") stop("fuzzy_crosspoints.py not found in package")
  reticulate::source_python(py)

  # 调用Python函数时传递用户手动输入的min_length, min_bin_size 和 seq_depth
  crosspoints_fuzzy(
    input_path = input_path, 
    output_path = output_path, 
    min_length = as.integer(min_length), 
    min_bin_size = as.integer(min_bin_size),
    seq_depth = as.numeric(seq_depth)
  )
}