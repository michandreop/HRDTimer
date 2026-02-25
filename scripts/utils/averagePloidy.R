averagePloidy <- function(bb) {
    c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
    sum(width(bb) * c * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}