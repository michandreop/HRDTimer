findMainCluster <- function(bb, min.dist=0.05){
    w <- which(bb$n.snv_mnv > 20 & !is.na(bb$time))
#   d <- dist(bb$time[w])
#   ci <- weighted.mean((bb$time.up - bb$time.lo)[w], width(bb)[w])
#   h <- hclust(d, method='average', members=bb$n.snv_mnv[w])
#   c <- cutree(h, h=ci)
#   ww <- c==which.max(table(c))
#   weighted.mean(bb$time[w][ww], 1/((bb$time.up - bb$time.lo + min.dist)[w][ww]), na.rm=TRUE)
    s <- seq(0,1,0.01)
    l2 <- pmin(bb$time.lo, bb$time - min.dist)[w]
    u2 <- pmax(bb$time.up, bb$time + min.dist)[w]
    l1 <- (l2 +  bb$time[w])/2
    u1 <- (u2+  bb$time[w])/2
    wd <- as.numeric(width(bb)[w])
    o <- sapply(s, function(i) sum(wd * ( (l2 <= i & u2 >=i) + (l1 <= i & u1 >= i))))
    s[which.max(o)]
}