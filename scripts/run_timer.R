 if (interactive()) {
    # default values for debugging mode
    # long sample ID

    aliquot_id <- '047f9e4d-86b5-4943-aef5-68199bf29e8c_p0.02'
    # purity
    purity <- as.numeric('0.639')
    data.type <- 'snvs'
    result_path <- "/home/dg204/park_dglodzik/TimeR_bb_cnv_test/"
    vcf_fn <- "/home/dg204/park_data/ICGC/SNV_indel_calls/final_consensus_12oct_passonly/snv_mnv/01658141-8398-4585-9f0f-8355dd9b0604.consensus.20160830.somatic.snv_mnv.vcf.gz"
    # file path to the copy number profile
    file <- "/home/dg204/park_data/ICGC/CNV_all_callers_included/01658141-8398-4585-9f0f-8355dd9b0604.consensus.20170119.somatic.cna.annotated.txt"
    genomeBuild <- 'hg37'
    tumor.id <- ''
     no_bootstraps <- 100
} else {
    args <- commandArgs(trailing = T)
    aliquot_id <- args[1]
    purity <-  as.numeric(args[2])
    data.type <- 'snvs'
    result_path <- args[3]
    vcf_fn <-  args[4]
    # file path to the copy number profile
    file <- args[5]
    genomeBuild <- args[6]
    tumor.id <- args[7]
    no_bootstraps <- 1000
}

data_fp <- paste0(result_path, "/data/")
if (!file.exists(data_fp)) {
  dir.create(data_fp, recursive = TRUE)
}
# will be used the create the RData object

vcf_fp <- paste0(result_path, "/vcfs/")
if (!file.exists(vcf_fp)) {
  dir.create(vcf_fp, recursive = TRUE)
}
# will be used the save the vcf object
plots_fp <- paste0(result_path, "/plots/")
if (!file.exists(plots_fp)) {
  dir.create(plots_fp, recursive = TRUE)
}
# for the Timer and copy number plot

# dependencies
suppressPackageStartupMessages({
    library(signature.tools.lib)
    library("extraDistr")
    library(TailRank)
    library("BSgenome.Hsapiens.UCSC.hg38")
    require(MutationTimeR)
})

source('utils/loadConsensusCNA.R')
source('utils/plotBB.R')
source('utils/plotSampleCustom.R')
source('utils/sampleVAFPlot.R')

# calculate parameters of beta-binomial distribution, given 
# mu : expected VAF (variant allele fraction)
# rho : overdispersion
# size : total number of observations (sequencing depth)
calculate_parameters <- function(mu, rho, size) {
  totalAB = size * (1 - rho) / rho
  alpha = mu * totalAB
  beta = (1 - mu) * totalAB
  return(list(alpha = alpha, beta = beta))
}


# Check version of the genome and build the catalogue
# this looks up mutation context for each mutation
if (genomeBuild=='hg38') {
    sample_muts <- vcfToSNVcatalogue(vcf_fn, genome.v="hg38")
} else if ((genomeBuild=='') || (genomeBuild=='hg37')) {
    sample_muts <- vcfToSNVcatalogue(vcf_fn)
}
sample_muts$muts$id <- paste0(sample_muts$muts$chroms, ':', sample_muts$muts$starts, '_', sample_muts$muts$wt, '/', sample_muts$muts$mt)
rownames(sample_muts$muts) <- sample_muts$muts$id 

# load the VCF file with variants
if (genomeBuild=='hg38') {
    vcf <- readVcf(vcf_fn, genome='hg38')
} else if ((genomeBuild=='') || (genomeBuild=='hg37')) {
    vcf <- readVcf(vcf_fn)
}
vcf@info$context <- sample_muts$muts[names(vcf),'context']
# prepare to add additional values to the VCF file
i = info(header(vcf))
newInfoHeader <- data.frame(Number=c(1),
                           Type=c('String'),
                           Description=c('mutation context'))
rownames(newInfoHeader) <- c('context')  
info(header(vcf)) <- rbind(i, newInfoHeader)

# add the allele counts
if (tumor.id!='') {
    counts_df <- data.frame(t_ref_count=sapply(geno(vcf)$AD[,tumor.id], head, 1), 
                            t_alt_count=sapply(geno(vcf)$AD[,tumor.id], tail, 1))
    i = info(header(vcf))
    newInfoHeader <- data.frame(Number=c(1,1),
                               Type=c('Integer', 'Integer'),
                               Description=c('tumor reference read count', 
                                             'tumor alternative read count'))
        
    rownames(newInfoHeader) <- c('t_ref_count','t_alt_count')  
    info(header(vcf)) <- rbind(i, newInfoHeader)
    info(vcf) <- cbind(info(vcf), counts_df)
}

# add the VAF field
if (! 'VAF' %in% colnames(info(vcf))) {
    # add variant allele fraction (VAF)
    i = info(header(vcf))
    newInfoHeader <- data.frame(Number=c(1),
                               Type=c('Float'),
                               Description=c('Variant Allele Fraction'))
        
    rownames(newInfoHeader) <- c('VAF')  
    info(header(vcf)) <- rbind(i, newInfoHeader)
    info(vcf)$VAF <- info(vcf)$t_alt_count/(info(vcf)$t_alt_count + info(vcf)$t_ref_count)
}



# if there are any variants
if (length(vcf)>0) {

    # load the copy number profile, depending on the format used
    if  (grepl('.consensus.20170119.somatic.cna.annotated.txt', file)) {
        # load and plot the copy number profile
        # tab-delimited file
        # chromosome, start, end,
        # total_cn, battenberg_nMaj2_A, battenberg_nMin2_A, battenberg_frac2_A, battenberg_nMaj1_A, battenberg_nMin1_A
        # battenberg_frac2_A
        # minor_cn, minor_cn 
        # battenberg_frac1_A, battenberg_frac2_A
        bb <- loadConsensusCNA(file, purity=purity)
        pdf(paste(plots_fp, aliquot_id, ".bb.pdf", sep=""), width=14, height=6)
        plotBB(bb)
        dev.off()   
    } else if (grepl('.ascat_ngs.summary.csv', file)) {
        tab <- read.table(file, header=FALSE, sep=',')
        gr <- GRanges(tab$V2, IRanges(tab$V3, tab$V4), strand="*", clonal_frequency=purity)
        seqlevels(gr) <- c(1:22, "X","Y")
        gr$major_cn <- tab$V7-tab$V8
        gr$minor_cn <- tab$V8
        gr$total_cn <- tab$V7
        bb <- gr
    } else if (grepl('.purple.segment.tsv', file) || grepl('.tsv', file)) {
        tab <- read.table(file, header=TRUE, sep='\t')
        gr <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", clonal_frequency=purity)
        seqlevels(gr) <- paste0('chr',c(1:22, "X","Y"))
        gr$major_cn <- round(tab$majorAlleleCopyNumber)
        gr$minor_cn <- round(tab$minorAlleleCopyNumber)
        gr$total_cn <- gr$major_cn + gr$minor_cn
        bb <- gr
    } else if (grepl('.RData', file)) {
        load(file)
    }
    
    # this is the main work of the TimeR package
    # 240611 - removed  isWgd=TRUE, as some of the samples we are analyzing do not have WGD, necessarily
    mt <- mutationTime(vcf, bb, n.boot=no_bootstraps)
    # add the inferred mutation time to the vcf file
    vcf <- addMutTime(vcf, mt$V)
    mcols(bb) <- cbind(mcols(bb),mt$T)

    # add power to detect mutations
    T <- vcf@info$MajCN + vcf@info$MinCN

    # estimates of mutation detection power (not used by HRDTimer directly
    rho=0.01 # set overdispersion
    # total coverage per mutation
    N_vec <- vcf@info$t_ref_count +  vcf@info$t_alt_count
    # preparing data frame for storing predicted mutation detection power
    power_df <- data.frame(N = N_vec, purity_adj = purity)
    # CNF will take value of purity, or fraction of clone with subclonal copy number change
    # fi : expected variant allele fraction of a mutation, given purity, total copy number (T) (also accounts for subclonal mutations through CNF)
    power_df$fi = (vcf@info$CNF * vcf@info$MutCN)/(purity * T + (1-purity)*2) 
    rowFunction <- function(row, rho) {
        params <- calculate_parameters(as.numeric(row['fi']), rho, as.numeric(row['N']))
        # probability that a mutation would be detected with 0, 1, or 2 reads, below what is detectable
        probability = dbb(0:2,  as.numeric(row['N']), params$alpha, params$beta)
        # probability that a mutation would be detcted
        pwr <- (1-sum(probability))
    }
    powr <- apply(power_df, 1, rowFunction,rho)
    vcf@info$fi <- power_df$fi
    vcf@info$powr <- powr
    # add info about new fields to the header
    newInfoHeader <- data.frame(Number=c(1,1),
                               Type=c('Float', 'Float'),
                               Description=c('expected allele fraction given inferred mutation CN', 
                                             'power to detect, given coverage at locus and expected VAF'))
    rownames(newInfoHeader) <- c('fi','powr')  
    i = info(header(vcf))
    info(header(vcf)) <- rbind(i, newInfoHeader)
    # the mutation detection power section is not directly used by HRDTimer
    # it was created to evaluate the possiblity that late mutations with lower VAF would be missed
    # but we found out that for most samples, the power to detect them is close to 100%
    
    # save the intermediate results
    save(mt, bb, vcf, file=paste(data_fp, aliquot_id, ".RData", sep=""))

    # make MutationTimeR plot as a sample-level summary
    if ((genomeBuild=='') || (genomeBuild=='hg37')) {
        pdf(paste(plots_fp, aliquot_id, ".pdf", sep=""))
        plotSample(vcf,bb)
        dev.off()
    } else if (genomeBuild=='hg38') {
        gr <- GRanges(
            seqnames = names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)),
            ranges = IRanges(rep(1,length(BSgenome.Hsapiens.UCSC.hg38)), 
                         end = seqlengths(BSgenome.Hsapiens.UCSC.hg38)),
            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
            )
        regions <- gr[seqnames(gr) %in% c(paste0('chr',1:22), 'chrX', 'chrY', 'chrM')]
        seqlevels(regions) <- c(paste0('chr',1:22), 'chrX', 'chrY', 'chrM')
        chrOffset <- c(0,cumsum(as.numeric(seqlengths(regions))[1:(length(seqlengths(regions))-1)]))
        names(chrOffset) <- names(seqlengths(regions))
        
        pdf(paste(result_path, aliquot_id, ".pdf", sep=""))
        plotSampleCustom(vcf,bb, regions=regions, chrOffset=chrOffset)
        dev.off()    
    }
    # end of the timeR plot

    # make a VAF plot
    sampleVAFPlot(p=purity, #p
                  vcf_df=as.data.frame(vcf@info), #vcf_df
                  writePdf=TRUE, # writePdf
                  plotFolder=plots_fp, # plotFolder
                  sample_id=aliquot_id)
    # end of VAF plot
    
    # write the vcf file with additional columns
    vcf_out_fn <- paste(vcf_fp, aliquot_id, ".vcf", sep="")
    writeVcf(vcf, vcf_out_fn)
    system(paste0('bgzip -c ', vcf_out_fn, ' > ', vcf_out_fn, '.gz'))
    system(paste0('tabix -p vcf ', vcf_out_fn, '.gz'))
    options(repr.plot.width=22, repr.plot.height=6)


    print('TimeR analyses complete')
}