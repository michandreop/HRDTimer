# this function concatenates clonal and subclonal segments
# keep track of copy number (clonal or subclonal in major_cn and minor_cn columns
# purity is the product of cancer cell fraction and clone frequency

loadConsensusCNA <- function(file, purity=1){
    if(grepl(".gz", file))
        file <- gzfile(file)
    tab <- read.table(file, header=TRUE, sep='\t')
    # non-missing values in the subclonal columns
    # defined wrt. tab
    subclonalIndex <- !is.na(tab$total_cn) & !is.na(tab$battenberg_nMaj2_A) & !is.na(tab$battenberg_nMin2_A) & !is.na(tab$battenberg_frac2_A) & (tab$battenberg_nMaj1_A == tab$major_cn & tab$battenberg_nMin1_A == tab$minor_cn | tab$battenberg_nMaj2_A == tab$major_cn & tab$battenberg_nMin2_A == tab$minor_cn)
    
    
    ix <- c(1:nrow(tab), which(subclonalIndex))
    gr <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", clonal_frequency=purity, tab[-3:-1])[ix]
    # gr have more rows than tabs
    # subclonal rows are being added to the clonal ones

       if(any(subclonalIndex)){
        # operations on the clonal rows
        gr$clonal_frequency[which(subclonalIndex)] <- tab$battenberg_frac1_A[subclonalIndex] * purity
        gr$major_cn[which(subclonalIndex)] <- tab$battenberg_nMaj1_A[subclonalIndex]
        gr$minor_cn[which(subclonalIndex)] <- tab$battenberg_nMin1_A[subclonalIndex]
        gr$total_cn[which(subclonalIndex)] <- tab$battenberg_nMaj1_A[subclonalIndex] + tab$battenberg_nMin1_A[subclonalIndex]

        # operations on the subclonal rows
        gr$clonal_frequency[-(1:nrow(tab))] <- tab$battenberg_frac2_A[subclonalIndex] * purity
        gr$total_cn[-(1:nrow(tab))] <- tab$battenberg_nMaj2_A[subclonalIndex] + tab$battenberg_nMin2_A[subclonalIndex]
        gr$major_cn[-(1:nrow(tab))] <- tab$battenberg_nMaj2_A[subclonalIndex]
        gr$minor_cn[-(1:nrow(tab))] <- tab$battenberg_nMin2_A[subclonalIndex]
    }
    seqlevels(gr) <- c(1:22, "X","Y")
    sort(gr)
}