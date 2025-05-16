    if (interactive()) {
    # default values for debugging mode
    # long sample ID

    aliquot_id <- 'fc8130df-2628-304a-e040-11ac0d485dfa'
    rdata_fn <- '~/park_dglodzik/TimeR_bb_all_v2//data/fc8130df-2628-304a-e040-11ac0d485dfa.RData'
    plot_folder <- '~/park_dglodzik/TimeR_bb_all_v2/hrdtimer_plots/'
    result_folder <- '~/park_dglodzik/TimeR_bb_all_v2/hrdtimer_results/'
    purity <- 0.461
     
        
    #aliquot_id <- '54e78de0-d357-4125-a904-ab35e461572b'
    # short sample ID
    #rdata_fn <- '/home/dg204/park_dglodzik/TimeR_bb_all_v2//data/54e78de0-d357-4125-a904-ab35e461572b.RData'
    #plot_folder <- '~/park_dglodzik/TimeR_bb_all_v2/hrdtimer_plots/'
    #result_folder <- '~/park_dglodzik/TimeR_bb_all_v2/hrdtimer_results/'
    #purity <- 0.461

    #aliquot_id <- '0d0793c1-df1b-4db1-ba36-adcb960cc0f5'
    #rdata_fn <- '/home/dg204/park_dglodzik/TimeR_bb_all_v2//data/0d0793c1-df1b-4db1-ba36-adcb960cc0f5.RData'
    #plot_folder <- '/home/dg204/park_dglodzik/TimeR_bb_all_v2//hrdtimer_plots/'
    #result_folder <- '/home/dg204/park_dglodzik/TimeR_bb_all_v2//hrdtimer_results/'
    #purity <- 0.74    
    
    #aliquot_id <- 'Patient101'
    # short sample ID
    #rdata_fn <- '/home/dg204/park_dglodzik/TimeR_INFORM/data/Patient101.RData'
    #plot_folder <- '/home/dg204/park_dglodzik/TimeR_INFORM/hrdtimer_plots/'
    #result_folder <- '/home/dg204/park_dglodzik/TimeR_INFORM/'
    #purity <- 0.68
    
    #aliquot_id <- 'AOCS-092-3-3:vs:AOCS-092-5-0'
    #rdata_fn <- "/home/dg204/park_dglodzik/TimeR_ICGC/data/AOCS-092-3-3:vs:AOCS-092-5-0.RData"
    #plot_folder <- "/home/dg204/park_dglodzik/TimeR_ICGC/hrdtimer_plots/"
    #result_folder <- "/home/dg204/park_dglodzik/TimeR_ICGC/hrdtimer_results/"
    #purity <- 0.59

    #aliquot_id <- 'PD31031a'
    #rdata_fn <- "/home/dg204/park_dglodzik/TimeR_SCANB/data/PD31031a.RData"
    #plot_folder <- "//home/dg204/park_dglodzik/TimeR_SCANB/hrdtimer_plots/"
    #result_folder <- "/home/dg204/park_dglodzik/TimeR_SCANB/hrdtimer_results/"
    #purity <- 0.31    

    nboot <- 200
    
} else {
    args <- commandArgs(trailing = T)
    aliquot_id <- args[1]
    rdata_fn <- args[2]
    plot_folder <- args[3]
    result_folder <- args[4]
    purity <-  as.numeric(args[5])
    nboot <- 20
}

print(paste('aliquot ID:', aliquot_id))
print(paste('RData timer:', rdata_fn))
print(paste('plot folder:', plot_folder))
print(paste('result folder:', result_folder))
print(paste('purity:', purity))

onlyClonal <- TRUE
selSig <- 'GEL-Breast_common_SBS3'

library(signature.tools.lib)
library(MutationTimeR)
library(Hmisc)
source('~/projects/rsignatures/src/utils/findMainCluster.R')
source('~/projects/brca_timing/src/utils/calculate_early_late_time.R')
source('~/projects/brca_timing/src/utils/calculate_early_late_time_prob.R')
source('~/projects/rsignatures/src/utils/plotBB.R')
source('~/projects/rsignatures/src/utils/plotDetailedTiming.R')
source('~/projects/brca_timing/src/utils/averageHom.R')
source('~/projects/brca_timing/src/utils/averagePloidy.R')
isDeamination <- function(vcf) grepl('.\\[C>T\\]G',info(vcf)$context)
source('~/projects/brca_timing/src/utils/cnTiming.R')
a <- function() {
    display(head(l$T[,c('type', 'time', 'time.lo' ,'time.up')]))
}
# 

sample.summary <- list()

#prop.clonal
#median.subclonal
#no.sbs.snps
# q1
# q5
# averageHom
# averagePloidy
# prop_genome_with_gains
# prop.muts.main.cn.cluster
# muts_sbs1_subclonal
# muts_sbs1_total
# G_deam
# muts_deam

if (file.exists(rdata_fn)) {
    
    load(rdata_fn)

    sample.summary[['isWGD']] <- classWgd(bb)
    
    #classWgd <- function(cn) .classWgd(getPloidy(cn), getHomozygousity(cn))
    
    # pre-process bb - the copy number file
    if (('battenberg_frac2_A' %in% colnames(bb)) && onlyClonal) {
        bb <- bb[is.na(bb$battenberg_frac2_A)]
    }
    bb$width <- end(bb) - start(bb)
    sample.summary[['prop.clonal']] <- sum(bb$width[!is.na(bb$battenberg_frac2_A)])/3e9
    sample.summary[['median.subclonal']] <- median(pmin(bb$battenberg_frac1_A, bb$battenberg_frac2_A), na.rm=TRUE)
    sample.summary[['genome.majorCN.bp']] <- sum(width(reduce(bb[!is.na(bb$major_cn) & (bb$major_cn==2)]))) 
    
    no.sbs.snps <- length(vcf)    

    # the mutation data frame
    info.df <- as.data.frame(info(vcf))
    info.df$chromosome <- as.character(seqnames(rowRanges(vcf)[rownames(info.df),]))
    # mutation cn status given total copy number
    info.df$mut.class <- paste(info.df$CLS, round(info.df$MutCN,1), '/', info.df$MajCN   )
    # allele fraction depending on on mutation copy number status

    info.df$tCN <- info.df$MajCN + info.df$MinCN

    # checking if the probability distribution over mutation copy number is well behaved
    q1 <- mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE)
    sample.summary[['q1']] <- q1
    q5 <- mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE)
    sample.summary[['q5']] <- q5
    n <- length(vcf)

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
    
    # proportion of genome with gains +1
    sample.summary[['prop_genome_with_gains']] <- round(sum(sum(coverage(bb[!is.na(bb$time.star) & bb$time.star=='***'])))/(3e9),2)

    # mutation where major copy number is 2 (minor copy number could be 1 or 2)
    info.df.majCN <- subset(info.df, MajCN==2 & (MutCN-round(MutCN,0))==0)

    # for internal purposes
    vaf_na <- is.na(info.df.majCN$VAF)
    if (sum(vaf_na)>0) {
        info.df.majCN$VAF[vaf_na] <-   info.df.majCN$t_alt_count[vaf_na] / rowSums(info.df.majCN[vaf_na,c('t_alt_count', 't_ref_count')])
    }
    info.df.majCN$u <- info.df.majCN$VAF * purity^(-1) * (purity * (info.df.majCN$MajCN + info.df.majCN$MinCN) + (1-purity) * 2)
    info.df.majCN$ccf <- info.df.majCN$u/ info.df.majCN$MutCN
    pdf(paste0(plot_folder, aliquot_id, '_ccf.pdf'), width=18, height=7)
    hist(info.df.majCN$ccf, breaks=100, border=NA, xlab='CCF (cancer cell fraction)')
    dev.off()
    
    
    if (length(unique(info.df.majCN$max_sig))>1) {

        # 1. Find segments compatible with WGD
        min.dist <- 0.05
        m <- findMainCluster(bb)
        print(paste('timing of the main cluster',m))
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

        i = info(header(vcf))
        newInfoHeader <- data.frame(Number=c(1),
                           Type=c('String'),
                           Description=c('whether in clean WGD region'
                                        ))
    
        rownames(newInfoHeader) <- c('inWGDregion')  
        info(header(vcf)) <- rbind(i, newInfoHeader)
        info(vcf)$inWGDregion <- as.character(vcf %over% bb[o])
        vcf_out_fn <- paste(result_folder, aliquot_id, ".vcf", sep="")
        writeVcf(vcf, vcf_out_fn)

        
         # 2. Find deaminations in compatible segments
        # removed the deamination part
        # deaminations in the major segment cluster, overlapping with 1 copy number segments, major copy number is 2
        # if the one below is selected, 
        w <- which(info(vcf)$MajCN==2 & sapply(info(vcf)$CNID, length)==1 &  vcf %over% bb[o] & isDeamination(vcf))

        
        
        sample.summary[['muts_sbs1_subclonal']] <- sum(info(vcf[w])$CLS=='subclonal')
        sample.summary[['muts_sbs1_total']] <- length(vcf[w])

        # if there is sufficient number of informative mutations
        if (length(w)>100) {

            # describe the signature
            result_path <- "~/projects/rsignatures/data/processed/TimeR_bb/"

            # make a signature plot for the molecular clock mutations
            #vcf_out_fn <- paste(result_folder, "/vcfs_sbs1/", aliquot_id, ".vcf", sep="")

            # over 2500 mutations, deaminations and not
            # mutations with major CN of 2, in the main cluster
            v <- vcf[w]
            print(paste(length(v), 'mutations with major CN of 2, in the main cluster'))
            if(nrow(v)<=90) return(NULL) # At least 100 SNVs
            if ('chr1' %in% seqlevels(v)) {
                seqnames(rowRanges(v)) <- factor(paste0('chr',3-info(v)$MinCN), levels=seqlevels(v))
            } else {
                seqnames(rowRanges(v)) <- factor(3-info(v)$MinCN, levels=seqlevels(v))
            }

            v <- sort(v)
            v_info <- info(v)
            v_info_df <- as.data.frame(v_info)
            v_info_df$state <- paste0(v_info_df$MajCN, '_', v_info_df$MinCN)
            v_info_df$isDeamination <-  isDeamination(v)
            
            all_vcf_w_info_df <- as.data.frame(info(vcf))
            all_vcf_w_info_df$isDeamination <-  isDeamination(vcf)
            
            all_vcf_w_info_df_deaminations <- subset(all_vcf_w_info_df, isDeamination)
            # normalize the mutation counts wrt ploidy, by looking up their copy number status
            sample.summary[['G_deam']] <- nrow(all_vcf_w_info_df_deaminations) / 
                            sum(all_vcf_w_info_df_deaminations$pMutCN / 
                            (all_vcf_w_info_df_deaminations$MajCN + all_vcf_w_info_df_deaminations$MinCN ), na.rm=TRUE)
            
            sample.summary[['muts_deam']] <- nrow(all_vcf_w_info_df_deaminations)
            
            # 3. Merged CN segments, by copy number
            # the approach is to merge consistent copy number segments
            # estimate the time of whole-genome duplication
            if ('chr1' %in% seqlevels(v)) {
                b <- GRanges(paste0('chr',1:3), IRanges(rep(1,3),rep(max(end(v)),3)), copy_number=4:2, major_cn=2, minor_cn=2:0, clonal_frequency=as.numeric(purity))
            } else {
                b <- GRanges(1:3, IRanges(rep(1,3),rep(max(end(v)),3)), copy_number=4:2, major_cn=2, minor_cn=2:0, clonal_frequency=as.numeric(purity))
            }
            l <- mutationTime(v, b, n.boot=nboot, isWgd=TRUE, rho=0.01, xmin=3)
            # the answer

            lt.df <- as.data.frame(l$T)
            rownames(lt.df) <- lt.df$type
            sample.summary[['timeR_2_2']] <- lt.df['Bi-allelic Gain (WGD)', 'time']
            sample.summary[['timeR_2_2_up']] <- lt.df['Bi-allelic Gain (WGD)', 'time.up']
            sample.summary[['timeR_2_2_lo']] <- lt.df['Bi-allelic Gain (WGD)', 'time.lo']
            sample.summary[['timeR_2_1']] <- lt.df['Mono-allelic Gain', 'time']
            sample.summary[['timeR_2_1_up']] <- lt.df['Mono-allelic Gain', 'time.up']
            sample.summary[['timeR_2_1_lo']] <- lt.df['Mono-allelic Gain', 'time.lo']
            sample.summary[['timeR_2_0']] <-  lt.df['CN-LOH', 'time']
            sample.summary[['timeR_2_0_up']] <- lt.df['CN-LOH', 'time.up']
            sample.summary[['timeR_2_0_lo']] <- lt.df['CN-LOH', 'time.lo']        
            sample.summary[['timeR']] <- weighted.mean(lt.df$time, lt.df$n.snv_mnv/sum(lt.df$n.snv_mnv))

            # molecular clock mutations only
            # Doga's method
            doga_timing <- calculate_early_late_time(v_info_df)
            doga_timing$ratio <- doga_timing$early/(doga_timing$early + doga_timing$late)
            rownames(doga_timing) <- doga_timing$state

            sample.summary[['doga_timing_wgd_2_2']] <- doga_timing['2_2', 'ratio']
            sample.summary[['doga_timing_wgd_2_1']] <- doga_timing['2_1', 'ratio'] 
            sample.summary[['doga_timing_wgd_2_0']] <- doga_timing['2_0', 'ratio']
            sample.summary[['doga_timing_wgd']] <- weighted.mean(doga_timing$ratio, doga_timing$nmuts/sum(doga_timing$nmuts))

            doga_timing_prob <- calculate_early_late_time_prob(v_info_df)
            doga_timing_prob$ratio <- doga_timing_prob$early/(doga_timing_prob$early + doga_timing_prob$late)
            rownames(doga_timing_prob) <- doga_timing_prob$state
            sample.summary[['doga_prob_timing_wgd']] <- weighted.mean(doga_timing_prob$ratio, doga_timing_prob$nmuts/sum(doga_timing_prob$nmuts))            
            
            w_sbs3 <- which(info(vcf)$MajCN==2 & sapply(info(vcf)$CNID, length)==1 &  vcf %over% bb[o] & info(vcf)$max_sig==selSig)
            # if there is a sufficient number of SBS3 alterations

            sample.summary[['muts_sbs3_early_chroms']] <- length(unique(seqnames(vcf[!is.na(info(vcf)$CLS) & info(vcf)$CLS=='clonal [early]',])))
            sample.summary[['muts_sbs3_subclonal']] <- sum(info(vcf[w_sbs3])$CLS=='subclonal')
            sample.summary[['muts_sbs3_total']] <- length(vcf[w_sbs3])

            # if there are at least 100 SBS3 mutations
            if (length(w_sbs3)>100) {
                print(paste(length(w_sbs3), ' SBS3 mutations ')) 
                
                all_vcf_w_info_df_sbs3 <- subset(all_vcf_w_info_df, max_sig=='SBS3')
                G_sbs3 <- nrow(all_vcf_w_info_df_sbs3) / sum(all_vcf_w_info_df_sbs3$pMutCN / (all_vcf_w_info_df_sbs3$MajCN + all_vcf_w_info_df_sbs3$MinCN ), na.rm=TRUE)
                # mutational signature exposure, corrected by average muation copy number
                sample.summary[['G_sbs3']] <- G_sbs3
                
                doga_timing_sbs3 <- calculate_early_late_time( as.data.frame(info(vcf[w_sbs3])))
                doga_timing_sbs3$ratio <- doga_timing_sbs3$early/(doga_timing_sbs3$early + doga_timing_sbs3$late)
                rownames(doga_timing_sbs3) <- doga_timing_sbs3$state

                sample.summary[['muts_sbs1_2_2_early']] <- doga_timing['2_2','early']
                sample.summary[['muts_sbs1_2_2_late']] <- doga_timing['2_2','late']
                sample.summary[['muts_sbs1_2_1_early']] <- doga_timing['2_1','early']
                sample.summary[['muts_sbs1_2_1_late']] <- doga_timing['2_1','late']
                sample.summary[['muts_sbs1_2_0_early']] <- doga_timing['2_0','early']
                sample.summary[['muts_sbs1_2_0_late']] <- doga_timing['2_0','late']

                sample.summary[['muts_sbs3_2_2_early']] <- doga_timing_sbs3['2_2','early']
                sample.summary[['muts_sbs3_2_2_late']] <- doga_timing_sbs3['2_2','late']
                sample.summary[['muts_sbs3_2_1_early']] <- doga_timing_sbs3['2_1','early']
                sample.summary[['muts_sbs3_2_1_late']] <- doga_timing_sbs3['2_1','late']
                sample.summary[['muts_sbs3_2_0_early']] <- doga_timing_sbs3['2_0','early']
                sample.summary[['muts_sbs3_2_0_late']] <- doga_timing_sbs3['2_0','late']            


                c_2_0_sbs3_sbs1 <- doga_timing['2_0','late']/ doga_timing_sbs3['2_0','late']
                # confidence intervals for binomial proportions
                c_2_0_conf <- as.data.frame(binconf(x=doga_timing['2_0','late'], n=doga_timing_sbs3['2_0','late']))

                t_hrd_2_0 <- (doga_timing['2_0','early'] - doga_timing_sbs3['2_0','early'] * c_2_0_sbs3_sbs1) / (doga_timing['2_0','early'] + doga_timing['2_0','late']) 
                sample.summary[['timing_hrd_2_0']] <- t_hrd_2_0
                t_hrd_2_0_lower <- (doga_timing['2_0','early'] - doga_timing_sbs3['2_0','early'] * c_2_0_conf$Lower) / (doga_timing['2_0','early'] + doga_timing['2_0','late']) 
                t_hrd_2_0_upper <- (doga_timing['2_0','early'] - doga_timing_sbs3['2_0','early'] * c_2_0_conf$Upper) / (doga_timing['2_0','early'] + doga_timing['2_0','late']) 
                sample.summary[['timing_hrd_2_0_lower']] <- t_hrd_2_0_lower
                sample.summary[['timing_hrd_2_0_upper']] <- t_hrd_2_0_lower
                t_wgd_2_0 <-  doga_timing['2_0','early'] / (doga_timing['2_0','early'] + doga_timing['2_0','late'])

                
                c_2_1_sbs3_sbs1 <- doga_timing['2_1','late'] / doga_timing_sbs3['2_1','late']
                c_2_1_conf <- as.data.frame(binconf(x=doga_timing['2_1','late'], n=doga_timing_sbs3['2_1','late']))
                t_hrd_2_1 <- (doga_timing['2_1','early'] - doga_timing_sbs3['2_1','early'] * c_2_1_sbs3_sbs1) / (doga_timing['2_1','early'] + doga_timing['2_1','late'])
                t_hrd_2_1_lower <- (doga_timing['2_1','early'] - doga_timing_sbs3['2_1','early'] * c_2_1_conf$Lower) / (doga_timing['2_1','early'] + doga_timing['2_1','late'])
                t_hrd_2_1_upper <- (doga_timing['2_1','early'] - doga_timing_sbs3['2_1','early'] * c_2_1_conf$Upper) / (doga_timing['2_1','early'] + doga_timing['2_1','late'])
                sample.summary[['timing_hrd_2_1']] <- t_hrd_2_1
                sample.summary[['timing_hrd_2_1_lower']] <- t_hrd_2_1_lower
                sample.summary[['timing_hrd_2_1_upper']] <- t_hrd_2_1_upper
                t_wgd_2_1 <-  doga_timing['2_1','early'] / (doga_timing['2_1','early'] + doga_timing['2_1','late'])
                
                c_2_2_sbs3_sbs1 <- doga_timing['2_2','late'] / doga_timing_sbs3['2_2','late']
                c_2_2_conf <- as.data.frame(binconf(x=doga_timing['2_2','late'], n=doga_timing_sbs3['2_2','late']))
                t_hrd_2_2 <- (doga_timing['2_2','early'] - doga_timing_sbs3['2_2','early'] * c_2_2_sbs3_sbs1) / (doga_timing['2_2','early'] + doga_timing['2_2','late'])
                t_hrd_2_2_lower <- (doga_timing['2_2','early'] - doga_timing_sbs3['2_2','early'] * c_2_2_sbs3_sbs1) / (doga_timing['2_2','early'] + doga_timing['2_2','late'])
                t_hrd_2_2_upper <- (doga_timing['2_2','early'] - doga_timing_sbs3['2_2','early'] * c_2_2_sbs3_sbs1) / (doga_timing['2_2','early'] + doga_timing['2_2','late'])
                sample.summary[['timing_hrd_2_2']] <- t_hrd_2_2
                sample.summary[['timing_hrd_2_2_lower']] <- t_hrd_2_2_lower
                sample.summary[['timing_hrd_2_2_upper']] <- t_hrd_2_2_upper
                t_wgd_2_2 <-  doga_timing['2_2','early'] / (doga_timing['2_2','early'] + doga_timing['2_2','late'])
                
                sample.summary[['c_2_0_sbs3_sbs1']] <- c_2_0_sbs3_sbs1
                sample.summary[['c_2_1_sbs3_sbs1']] <- c_2_1_sbs3_sbs1
                sample.summary[['c_2_2_sbs3_sbs1']] <- c_2_2_sbs3_sbs1
                
                # now also set up HRD
                if ((!is.na(lt.df['Bi-allelic Gain (WGD)','time'])) & ((lt.df['Bi-allelic Gain (WGD)','time.up'] - lt.df['Bi-allelic Gain (WGD)','time.lo']) > 0.5)) {
                    t_hrd_2_2 <- NA
                    t_wgd_2_2 <- NA
                }
                if ((!is.na(lt.df['Mono-allelic Gain','time'])) & ((lt.df['Mono-allelic Gain','time.up'] - lt.df['Mono-allelic Gain','time.lo']) > 0.5)) {
                    t_hrd_2_1 <- NA
                    t_wgd_2_1 <- NA
                }
                if ((!is.na(lt.df['CN-LOH','time'])) & ((lt.df['CN-LOH','time.up'] - lt.df['CN-LOH','time.lo']) > 0.5)) {
                    t_hrd_2_0 <- NA
                    t_wgd_2_0 <- NA
                }                    

                sample.summary[['timing_hrd']] <- weighted.mean(c(t_hrd_2_2, t_hrd_2_1, t_hrd_2_0),
                                                                     c(lt.df['Bi-allelic Gain (WGD)', 'n.snv_mnv'],
                                                                        lt.df['Mono-allelic Gain', 'n.snv_mnv'], 
                                                                       lt.df['CN-LOH', 'n.snv_mnv'])/sum(lt.df$n.snv_mnv), na.rm = TRUE)

                sample.summary[['timing_wgd_doga']] <- weighted.mean(c(t_wgd_2_2, t_wgd_2_1, t_wgd_2_0),
                                                                     c(lt.df['Bi-allelic Gain (WGD)', 'n.snv_mnv'],
                                                                        lt.df['Mono-allelic Gain', 'n.snv_mnv'], 
                                                                       lt.df['CN-LOH', 'n.snv_mnv'])/sum(lt.df$n.snv_mnv), na.rm = TRUE)   

                
            } # end if SBS3

        } # end if SBS1               
    } # if there are at least 2 signatures
    save(sample.summary, file=paste0(result_folder, aliquot_id, '.RData'))
} # if file exists