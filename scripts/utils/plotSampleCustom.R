 plotSampleCustom <- function (vcf, cn, sv = NULL, title = "", regions = NULL, ylim.cn = c(0, 
    5), layout.height = c(4, 1.2, 3.5), y.sv = ylim.cn[2] - 1, 
    chrOffset= NULL) 
{
    #if (is.null(regions)) 
    #    regions <- refLengths[1:24]
     
    p <- par()
    layout(matrix(1:3, ncol = 1), height = layout.height)
    par(mar = c(0.5, 3, 0.5, 0.5), mgp = c(2, 0.25, 0), bty = "L", 
        las = 2, tcl = -0.25, cex = 1)
     
    xlim = c(min(chrOffset[as.character(seqnames(regions))] + 
        start(regions)), max(chrOffset[as.character(seqnames(regions))] + 
        end(regions)))
    
     bbb <- cn[cn %over% regions]


    plotVcfCustom(vcf[vcf %over% regions], bbb, legend = FALSE, col.grid = "white", 
        xaxt = FALSE, cex = 0.33, xlim = xlim)
     
    mtext(line = -1, side = 3, title, las = 1)
    plotBBCustom(bbb, 
                 ylim = ylim.cn, 
                 legend = FALSE, 
                 type = "bar", 
                 col.grid = "white", 
                 col = c("lightgrey", "darkgrey"), 
                 xaxt = FALSE, 
                 xlim = xlim,
                 refLengths=regions
        )
     
    par(mar = c(3, 3, 0.5, 0.5))
    plotTimingCustom(bbb, xlim = xlim, legend = FALSE, col.grid = NA, chrOffset=chrOffset)
    

    
    #if (any(!is.na(cn$time))) {
    #    y0 <- seq(0.005, 0.995, 0.01)
    #    s <- .histBeta(cn)
    #    g <- colorRampPalette(RColorBrewer::brewer.pal(4, "Set1")[c(3, 
    #        2, 4)])(100)
    #    segments(x0 = chrOffset["MT"], y0 = y0, x1 = chrOffset["MT"] + 
    #        s/max(s) * 1e+08, col = g, lend = 3)
    #    getMode <- function(s) {
    #        if (all(is.na(s))) 
    #            return(NA)
    #        w <- which.max(s)
    #        if (w %in% c(1, length(s))) {
    #            m <- which(c(0, diff(s)) > 0 & c(diff(s), 0) < 
    #              0)
    #            if (length(m) == 0) 
    #              return(w)
    #            m <- m[which.max(s[m])]
    #            return(if (s[w] > 2 * s[m]) w else m)
    #        }
    #        else return(w)
    #    }
    #    abline(h = y0[getMode(s)], lty = 5)
    #    if ("time.2nd" %in% colnames(mcols(cn))) 
    #        if (any(!is.na(cn$time.2nd))) {
    #            s2 <- .histBeta(cn, time = "time.2nd")
    #            segments(x0 = chrOffset["MT"] + s/max(s) * 1e+08, 
    #              y0 = y0, x1 = chrOffset["MT"] + s/max(s) * 
    #                1e+08 + s2/max(s) * 1e+08, col = paste0(g, 
    #                "44"), lend = 3)
    #            abline(h = y0[getMode(s2)], lty = 3)
    #        }
    #}

                       
    par(p[setdiff(names(p), c("cin", "cra", "csi", "cxy", "din", 
        "page"))])
}

                       
plotVcfCustom <- function(vcf, 
                          bb, 
                          col = RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)], 
                          ID = meta(header(vcf))[[1]]["ID",1], 
                          legend=TRUE, 
                          lty.grid=1, 
                          col.grid="grey", 
                          xaxt=TRUE, 
                          pch=16, 
                          pch.out=pch, 
                          cex=0.66, 
                          xlim=c(0,chrOffset["MT"])) {
    
	cls <- factor(paste(as.character(info(vcf)$CLS)), levels = c("clonal [early]","clonal [late]","clonal [NA]","subclonal" , "NA"))
	plot(NA,NA, xlab='', ylab="VAF", ylim=c(0,1), xlim=xlim, xaxt="n", cex=cex)
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	for(i in which(!sapply(bb$timing_param, is.null))) {
		s <- start(bb)[i]
		e <- end(bb)[i]
		x <- chrOffset[as.character(seqnames(bb)[i])]
		y <- bb$timing_param[[i]][,"f"]
		l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
		l[is.na(l)] <- 0
		if(any(is.na(c(s,e,x,y,l)))) next
		segments(s+x,y,e+x,y, lwd=l*4+.1)
	}
	points(start(vcf) + chrOffset[as.character(seqnames(vcf))], getAltCount(vcf)/getTumorDepth(vcf),col=col[cls],  pch=ifelse(info(vcf)$pMutCNTail < 0.025 | info(vcf)$pMutCNTail > 0.975, pch.out , pch),  cex=cex)				
	if(legend) legend("topleft", pch=19, col=col, legend=paste(as.numeric(table(cls)), levels(cls)), bg='white')
}

plotBBCustom <- function(bb, 
                         ylim=c(0,max(max(bb$total_cn, na.rm=TRUE))), 
                         col=RColorBrewer::brewer.pal(4,"Set2"), type=c("lines","bars"), 
                         legend=TRUE, 
                         lty.grid=1, 
                         col.grid="grey", 
                         xaxt=TRUE, xlim=c(min(chrOffset[as.character(seqnames(bb))]+start(bb)),max(chrOffset[as.character(seqnames(bb))]+end(bb))),
                        refLengths=refLengths){
    
	type <- match.arg(type)
	#s <- c(1:22, "X","Y")
    s <- seqnames(refLengths)
	l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
	names(l) <- s
	plot(NA,NA, ylab="Copy number",xlab="",xlim=xlim, ylim=ylim, xaxt="n")
	c <- cumsum(l)
	axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
	if(xaxt) mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
	#abline(v=c, lty=3)

   
	if(type=="lines"){
		x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
		x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
		lwd <- 5* bb$clonal_frequency / max(bb$clonal_frequency)
		segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=col[1], lwd=lwd)
		segments(x0=x0, bb$minor_cn -.125,x1, bb$minor_cn-.125, col=col[2], lwd=lwd)
		segments(x0=x0, bb$total_cn+.125,x1, bb$total_cn+.125, col=1, lwd=lwd)

	}else{
		ub <- unique(bb)
		f <- findOverlaps(ub,bb)

		m <- t(model.matrix( ~ 0 + factor(queryHits(f))))
		ub$major_cn <- m %*% mg14::na.zero(bb$major_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
		ub$minor_cn <- m %*% mg14::na.zero(bb$minor_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
		ub$total_cn <- ub$major_cn + ub$minor_cn
		ub$clonal_frequency <- max(bb$clonal_frequency)
		x0 <- start(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
		x1 <- end(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
		rect(x0,0,x1, ub$major_cn, col=col[2], lwd=NA)
		rect(x0,ub$major_cn,x1, ub$total_cn, col=col[1], lwd=NA)
		abline(h = 1:floor(ylim[2]), lty=lty.grid, col=col.grid)
	}
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	if(legend){
		if(type=="lines") legend("topleft", c("Total CN","Major CN","Minor CN"), col=c("black", col[1:2]), lty=1, lwd=2, bg='white')
		else legend("topleft", c("Major CN","Minor CN"), fill=col[1:2], bg='white')
	}
}

plotTimingCustom <- function(bb, 
                             time=mcols(bb), 
                             col=paste0(RColorBrewer::brewer.pal(5,"Set2")[c(3:5)],"88"), 
                             legend=TRUE, 
                             col.grid='grey', 
                             lty.grid=1, 
                             xlim=c(0,chrOffset["MT"]), 
                             plot=2,
                            chrOffset=chrOffset){
    
	plot(NA,NA, xlab='', ylab="Time [mutations]", ylim=c(0,1), xlim=xlim, xaxt="n")
	if(any(!is.na(bb$time)))
		tryCatch({
					bb <- bb[!is.na(bb$time)]
					s <- start(bb)
					e <- end(bb)
					x <- chrOffset[as.character(seqnames(bb))]
					y <- time[,"time"]
					rect(s+x,time[,"time.lo"],e+x,time[,"time.up"], border=NA, col=col[time[,"type"]], angle = ifelse(bb$time.star=="*" | is.na(bb$time.star),45,135), density=ifelse(bb$time.star == "***", -1, 72))
					segments(s+x,y,e+x,y)
					
					if("time.2nd" %in% colnames(time)){ 
						w <- !is.na(time[,"time.2nd"])
						if(sum(w) != 0 & plot==2){
							s <- start(bb)[w]
							e <- end(bb)[w]
							x <- chrOffset[as.character(seqnames(bb))][w]
							y <- time[w,"time.2nd"]
							rect(s+x,time[w,"time.2nd.lo"],e+x,time[w,"time.2nd.up"], border=NA, col=sub("88$","44",col)[as.numeric(time[w,"type"])], angle = ifelse(bb$time.star[w]=="*" | is.na(bb$time.star[w]),45,135), density=ifelse(bb$time.star[w] == "***", -1, 72))
							segments(s+x,y,e+x,y)
						}
					}
				}, error=function(x) warning(x))
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	s <- c(1:22, "X","Y")
	l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
	names(l) <- s
	c <- cumsum(l)
	axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
	mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	if(legend) legend("topleft", levels(time[,"type"]), fill=col, bg="white")
}
