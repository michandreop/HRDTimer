library(Rsamtools)
library(VariantAnnotation)

refFile = "/home/dg204/park_dglodzik/ref_cgpwgs/core_ref_GRCh37d5/genome.fa" #meta(header(v))["reference",]
refLengths <- scanFaIndex(file=refFile)
chrOffset <- cumsum(c(0,as.numeric(width(refLengths))))
names(chrOffset) <- c(seqlevels(refLengths), "NA")




plotBB <- function(bb, ylim=c(0,max(max(bb$total_cn, na.rm=TRUE))), 
                   col=RColorBrewer::brewer.pal(4,"Set2"), type=c("lines","bars"), legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, 
                   xlim=c(min(chrOffset[as.character(seqnames(bb))]+start(bb)),max(chrOffset[as.character(seqnames(bb))]+end(bb))), add=FALSE, regions=refLengths[1:24]){
    type <- match.arg(type)
    s <- c(1:22, "X","Y")
    l <- as.numeric(width(refLengths[seqnames(refLengths) %in% s]))
    names(l) <- s
    c <- cumsum(l)
    if(add==FALSE){
        plot(NA,NA, ylab="Copy number",xlab="",xlim=xlim, ylim=ylim, xaxt="n")
        axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
        if(length(regions) == 1)
            if(xaxt) axis(side=1, at=pretty(c(start(regions), end(regions)))+chrOffset[as.character(seqnames(regions))], labels=sitools::f2si(pretty(c(start(regions), end(regions)))), las=1)
        if(xaxt) mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
    }
    if(type=="lines"){
    x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
    x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
    lwd <- 5* bb$clonal_frequency / max(bb$clonal_frequency)
    segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=col[1], lwd=lwd)
    segments(x0=x0, bb$minor_cn -.125,x1, bb$minor_cn-.125, col=col[2], lwd=lwd)
    segments(x0=x0, bb$total_cn+.125,x1, bb$total_cn+.125, col=1, lwd=lwd)
#   cv <- coverage(bb)
#   cv <- cv[s[s%in%names(cv)]]
#   par(xpd=NA)
#   for(n in names(cv)){
#       cc <- cv[[n]]
#       segments(start(cc) + cumsum(l)[n] - l[n] ,-runValue(cc)/2,end(cc)+ cumsum(l)[n] - l[n], -runValue(cc)/2, col=4)
#   }
    }else{
        ub <- unique(bb)
        f <- findOverlaps(ub,bb)
        m <- t(model.matrix( ~ 0 + factor(queryHits(f))))
        ub$total_cn <- m %*% mg14::na.zero(bb$total_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
        ub$major_cn <- m %*% mg14::na.zero(bb$major_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
        ub$minor_cn <- m %*% mg14::na.zero(bb$minor_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
        ub$clonal_frequency <- max(bb$clonal_frequency)
        x0 <- start(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
        x1 <- end(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
        rect(x0,0,x1, ub$major_cn, col=col[2], lwd=NA)
        rect(x0,ub$major_cn,x1, ub$total_cn, col=col[1], lwd=NA)
        abline(h = 1:floor(ylim[2]), lty=lty.grid, col=col.grid)
    }
    abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
    #if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
    if(legend){
        if(type=="lines") legend("topleft", c("Total CN","Major CN","Minor CN"), col=c("black", col[1:2]), lty=1, lwd=2, bg='white')
        else legend("topleft", c("Major CN","Minor CN"), fill=col[1:2], bg='white')
    }
}