# looks for regions consistent with genome doubling
# annotates vcf files with with inWGDregion tag

if (interactive()) {
    # default values for debugging mode
    # long sample IDx

    aliquot_id <- '008aef39-0c97-48ce-9dfd-f12d67116c59_p0.5590'
    rdata_fn <- '/home/dg204/park_dglodzik/TimeR_bb_purity_exp//data/008aef39-0c97-48ce-9dfd-f12d67116c59_p0.5590.RData'
    plot_folder <- '/home/dg204/park_dglodzik/TimeR_bb_purity_exp//hrdtimer_plots/'
    result_folder <- '/home/dg204/park_dglodzik/TimeR_bb_purity_exp//hrdtimer_results/'
    purity <- 0.5590
    
} else {
    args <- commandArgs(trailing = T)
    aliquot_id <- args[1]
    rdata_fn <- args[2]
    plot_folder <- args[3]
    result_folder <- args[4]
    purity <-  as.numeric(args[5])
}

print(paste('aliquot ID:', aliquot_id))
print(paste('RData timer:', rdata_fn))
print(paste('plot folder:', plot_folder))
print(paste('result folder:', result_folder))
print(paste('purity:', purity))

# remove regions with subclonal copy number changes
onlyClonal <- TRUE

suppressPackageStartupMessages({
library(signature.tools.lib)
library(MutationTimeR)
library(Hmisc)
})

source('utils/findMainCluster.R')
source('utils/averageHom.R')
source('utils/averagePloidy.R')

sample.summary <- list()
# isWGD : whether the sample is predicted to have a genome doubling
# prop.clonal : proportion of genome (by copy number footprint) that is clonal
# median.subclonal : median subclonal fraction of sublconal CNV segments
# genome.majorCN.bp : genome footprint where major CN = 2 (amenable for timing analysis)
# q1 : tail probability of the VAF distribution (1%)
# q5 : tail probability of the VAF distribution (5%)
# averageHom : the proportion of the genome that is "homozygous" (where the minor allele count is zero)
# averagePloidy : average ploidy across copy number segments
# main.cluster.timing : WGD timing (based on all mutations, and not stratified by signatures)
# prop.muts.main.cn.cluster : proportion of the genome in the main timing cluster
# prop.muts.main.cn.cluster.majorCN : proportion of the genome in the main cluster, with majorCN=2
# meanMutCN : mean mutation copy number (multiplicity)

if (file.exists(rdata_fn)) {

    # loading the result of MutationTimeR 
    load(rdata_fn)

    sample.summary[['isWGD']] <- classWgd(bb)

    
    # pre-process bb - the copy number file
    if (('battenberg_frac2_A' %in% colnames(bb)) && onlyClonal) {
        # remove regions with subclonal copy number changes
        bb <- bb[is.na(bb$battenberg_frac2_A)]
    }
    # segment footprint (bp)
    bb$width <- end(bb) - start(bb)
    # proportion of genome (by copy number footprint) that is clonal
    sample.summary[['prop.clonal']] <- sum(bb$width[!is.na(bb$battenberg_frac2_A)])/3e9
    # median subclonal fraction of sublconal CNV segments
    sample.summary[['median.subclonal']] <- median(pmin(bb$battenberg_frac1_A, bb$battenberg_frac2_A), na.rm=TRUE)
    # genome footprint where major CN = 2 (amenable for timing analysis)
    sample.summary[['genome.majorCN.bp']] <- sum(width(reduce(bb[!is.na(bb$major_cn) & (bb$major_cn==2)]))) 

    # number of mutations
    no.sbs.snps <- length(vcf)    

    # the MUTATION data frame
    info.df <- as.data.frame(info(vcf))
    info.df$chromosome <- as.character(seqnames(rowRanges(vcf)[rownames(info.df),]))
    # mutation cn status given total copy number
    info.df$mut.class <- paste(info.df$CLS, round(info.df$MutCN,1), '/', info.df$MajCN   )
    # local tumor total copy number
    info.df$tCN <- info.df$MajCN + info.df$MinCN

    # checking if the probability distribution over mutation copy number is well behaved
    q1 <- mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE)
    sample.summary[['q1']] <- q1
    q5 <- mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE)
    sample.summary[['q5']] <- q5

    # annotating segments to those where the confidence interval is reasonable
    mt$T$is.certain <- !is.na(mt$T$time) & (mt$T$time.up - mt$T$time.lo)<0.5
    # this is actually a segment data frame
    mt.df <- as.data.frame(mt$T)

    certain.segments <- subset(mt$T, is.certain==TRUE )
    if (nrow(certain.segments)>0) {
        mt.df$certainty <- 'uncertain'
        mt.df$certainty[mt.df$is.certain] <- 'certain'
    }
    sample.summary[['averageHom']] <- averageHom(bb)
    sample.summary[['averagePloidy']] <- averagePloidy(bb)


    # mutation where major copy number is 2 (minor copy number could be 1 or 2)
    info.df.majCN <- subset(info.df, MajCN==2 & (MutCN-round(MutCN,0))==0)
    if (nrow(info.df.majCN)>0) {
        # for debugging purposes
        # calculate VAF
        vaf_na <- is.na(info.df.majCN$VAF)
        if (sum(vaf_na)>0) {
            info.df.majCN$VAF[vaf_na] <-   info.df.majCN$t_alt_count[vaf_na] / rowSums(info.df.majCN[vaf_na,c('t_alt_count', 't_ref_count')])
        }
        # calculate cancer cell fraction, based on inferred mutation multiplicity
        info.df.majCN$u <- info.df.majCN$VAF * purity^(-1) * (purity * (info.df.majCN$MajCN + info.df.majCN$MinCN) + (1-purity) * 2)
        info.df.majCN$ccf <- info.df.majCN$u/ info.df.majCN$MutCN
        # save the CCF profile
        pdf(paste0(plot_folder, aliquot_id, '_ccf.pdf'), width=18, height=7)
        hist(info.df.majCN$ccf, breaks=100, border=NA, xlab='CCF (cancer cell fraction)')
        dev.off()
    }


    # 1. Find segments compatible with WGD
    min.dist <- 0.05
    m <- findMainCluster(bb)
    print(paste('timing of the main cluster',m))
    sample.summary[['main.cluster.timing']] <- m
    # lower and upper boundaries for the 
    l <- pmin(bb$time.lo, bb$time - min.dist)
    u <- pmax(bb$time.up, bb$time + min.dist)
    # events overlapping with the main cluster
    o <- which(l <= m & u >= m)

    # proportion of the genome in the main timing cluster
    sample.summary[['prop.muts.main.cn.cluster']] <- round(sum(sum(coverage(bb[o])))/(3e9),2)

    # proportion of the genome in the main cluster, with majorCN=2
    bb_o <- bb[o]
    sample.summary[['prop.muts.main.cn.cluster.majorCN']] <- round(sum(sum(coverage(bb_o[bb_o$major_cn==2])))/(3e9),2)

    # add a new info field to the .vcf file
    i = info(header(vcf))
    newInfoHeader <- data.frame(Number=c(1),
                       Type=c('String'),
                       Description=c('whether in clean WGD region'
                                    ))

    rownames(newInfoHeader) <- c('inWGDregion')  
    info(header(vcf)) <- rbind(i, newInfoHeader)
    # this flags each mutation whether it had occurred in a 
    info(vcf)$inWGDregion <- as.character(vcf %over% bb[o])
    vcf_out_fn <- paste(result_folder, aliquot_id, ".vcf", sep="")
    writeVcf(vcf, vcf_out_fn)

    # mean mutation copy number
    sample.summary[['meanMutCN']] <-  mean(info(vcf)$MutCN, na.rm=TRUE)
    

    print(paste('Writing results to:', paste0(result_folder, aliquot_id, '.RData')))
    save(sample.summary, file=paste0(result_folder, aliquot_id, '.RData'))
} # if file exists