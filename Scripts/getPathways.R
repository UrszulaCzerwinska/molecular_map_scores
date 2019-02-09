gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

gmt_to_table <- function(gmt.list) {
  n.obs <- sapply(gmt.list, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(gmt.list, "[", i = seq.max))
  map.gmt = data.frame(t(mat), stringsAsFactors = TRUE)
  return(map.gmt)
}
  