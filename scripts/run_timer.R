 if (interactive()) {
    # default values for debugging mode
    # long sample ID
     
    aliquot_id <- 'PD31031a'
    # purity
    purity <- as.numeric('0.31')
    data.type <- 'snvs'
    result_path <- "~/park_dglodzik/TimeR_SCANB/"
    vcf_fn <- "/home/dg204/repos/indels_sigmatrix/data/vcf_files_subs/PD31031a_subs.vcf.gz"
    # file path to the copy number profile
    file <- "/home/dg204/park_dglodzik/data_repo/Lund/calls/Profiles/PD31031a.ascat_ngs.summary.csv"
    genomeBuild <- 'hg37'
    tumor.id <- ''
    no_bootstraps <- 10
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
vcf_fp <- paste0(result_path, "/vcfs/")
if (!file.exists(vcf_fp)) {
  dir.create(vcf_fp, recursive = TRUE)
}
plots_fp <- paste0(result_path, "/plots/")
if (!file.exists(plots_fp)) {
  dir.create(plots_fp, recursive = TRUE)
}
sig_fp <- paste0(result_path, "/sigs/")
if (!file.exists(data_fp)) {
  dir.create(data_fp, recursive = TRUE)
}


# dependencies
library(signature.tools.lib)
library("extraDistr")
library(TailRank)
source('utils/loadConsensusCNA.R')
source('utils/plotBB.R')
source('utils/plotSampleCustom.R')
library("BSgenome.Hsapiens.UCSC.hg38")

require(MutationTimeR)
#source('/home/dg204/repos/MutationTimeR/R/MutationTime.R')


calculate_parameters <- function(mu, rho, size) {
  totalAB = size * (1 - rho) / rho
  alpha = mu * totalAB
  beta = (1 - mu) * totalAB
  return(list(alpha = alpha, beta = beta))
}

print(aliquot_id)
print(purity)



# Check version of the genome and build the catalogie
if (genomeBuild=='hg38') {
    sample_muts <- vcfToSNVcatalogue(vcf_fn, genome.v="hg38")
} else if ((genomeBuild=='') || (genomeBuild=='hg37')) {
    sample_muts <- vcfToSNVcatalogue(vcf_fn)
}
sample_muts$muts$id <- paste0(sample_muts$muts$chroms, ':', sample_muts$muts$starts, '_', sample_muts$muts$wt, '/', sample_muts$muts$mt)
rownames(sample_muts$muts) <- sample_muts$muts$id 

options(repr.plot.width=22, repr.plot.height=6)
# plot mutational profile of a given sample
pdf(paste(plots_fp, aliquot_id, ".sigs.pdf", sep=""))
plotSubsSignatures(sample_muts$catalogue,  add_to_titles=aliquot_id)
dev.off()

# fit the known signatures
fit.result <- FitMS(sample_muts$catalogue,
                   organ='Breast')
plotFitMS(fit.result, outdir = paste0(result_path, aliquot_id, '/'))

allSignatures.m <- cbind(fit.result$commonSignatures,
                      fit.result$rareSignatures)

assessed.signatures <- colnames(fit.result$exposures)
assessed.signatures <- assessed.signatures[assessed.signatures!='unassigned']

allSignatures.m <- cbind(fit.result$commonSignatures,
                      fit.result$rareSignatures
                      )[,assessed.signatures]
fit.result$exposures.norm <- fit.result$exposures/sum(fit.result$exposures)
exposure.m <- matrix(fit.result$exposures.norm[,assessed.signatures],
                    nrow=nrow(allSignatures.m),
                    ncol=ncol(allSignatures.m), byrow=TRUE)
P_s_given_context <- allSignatures.m * exposure.m
P_s_given_context_norm<- P_s_given_context/rowSums(P_s_given_context)

max_sig_prob <- apply(P_s_given_context, 1, max)
max_sig <- assessed.signatures[max.col(P_s_given_context, "first")]
names(max_sig) <- rownames(allSignatures.m)

# load the VCF file with variants
if (genomeBuild=='hg38') {
    vcf <- readVcf(vcf_fn, genome='hg38')
} else if ((genomeBuild=='') || (genomeBuild=='hg37')) {
    vcf <- readVcf(vcf_fn)
}
vcf@info$context <- sample_muts$muts[names(vcf),'context']
vcf@info$max_sig <- max_sig[vcf@info$context ]

if ('GEL-Breast_common_SBS3' %in% assessed.signatures) {
    vcf@info$prob_SBS3 <- P_s_given_context_norm[vcf@info$context,'GEL-Breast_common_SBS3']  
} else {
    vcf@info$prob_SBS3 <- NA
}

if ('GEL-Breast_common_SBS1' %in% assessed.signatures) {
    vcf@info$prob_SBS1 <- P_s_given_context_norm[vcf@info$context,'GEL-Breast_common_SBS1']  
} else {
    vcf@info$prob_SBS1 <- NA
}



i = info(header(vcf))
newInfoHeader <- data.frame(Number=c(1,1,1,1),
                           Type=c('String', 'String', 'Float', 'Float'),
                           Description=c('mutation context',
                                         'most likely signature',
                                         'probability of SBS3',
                                         'probability of SBS1'
                                        ))
    
rownames(newInfoHeader) <- c('context', 'max_sig', 'prob_SBS3', 'prob_SBS1')  
info(header(vcf)) <- rbind(i, newInfoHeader)

if (tumor.id!='') {
        # add the allele counts
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

if (! 'VAF' %in% colnames(info(vcf))) {

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
    }
    
    # this is the main work of the TimeR package
    # 240611 - removed  isWgd=TRUE, as some of the samples we are analyzing do not have WGD, necessarily
    mt <- mutationTime(vcf, bb, n.boot=no_bootstraps)
    vcf <- addMutTime(vcf, mt$V)
    mcols(bb) <- cbind(mcols(bb),mt$T)

    # add power to dectect mutations
    T <- vcf@info$MajCN + vcf@info$MinCN
    
    
    rho=0.01 
    N_vec <- vcf@info$t_ref_count +  vcf@info$t_alt_count
    power_df <- data.frame(N = N_vec, purity_adj = purity)
    # CNF will take value of purity, or half of purity for subclonal mutations
    power_df$fi = (vcf@info$CNF * vcf@info$MutCN)/(purity * T + (1-purity)*2) 
    rowFunction <- function(row, rho) {
        params <- calculate_parameters(as.numeric(row['fi']), rho, as.numeric(row['N']))
        probability = dbb(0:2,  as.numeric(row['N']), params$alpha, params$beta)
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

    
    # save the intermediate results
    save(mt, bb, vcf, file=paste(data_fp, aliquot_id, ".RData", sep=""))

    # make timeR plot
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
    
    # write the vcf file with additional columns
    vcf_out_fn <- paste(vcf_fp, aliquot_id, ".vcf", sep="")
    writeVcf(vcf, vcf_out_fn)
    system(paste0('bgzip -c ', vcf_out_fn, ' > ', vcf_out_fn, '.gz'))
    system(paste0('tabix -p vcf ', vcf_out_fn, '.gz'))
    options(repr.plot.width=22, repr.plot.height=6)


    print('TimeR analyses complete')
}