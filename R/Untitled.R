
m <- 3000
simulate.p.values <- function(
    i, delta = 2, fraction = 0.1,ss = 4){
  control <- rnorm(ss,0,1)
  treatment <- rnorm(ss,0,1)
  if (runif(1) < fraction)
    treatment <- treatment + delta
  return(t.test(treatment,control)$p.value)
}
pvals00 <- sapply(1:m, simulate.p.values,
                  delta = 0,fraction = 0 )
pvals24 <- sapply(1:m, simulate.p.values,
                  delta = 2, fraction = 0.4 )

par(mfrow=c(1,3))
hist(pvals00, breaks=20, ylim=c(0,m/3), main="A")
abline(h=m/20, col=3)
abline(v = 0.05, col=2)
text(x = c(0.02,0.02,0.2,0.2),y=c(20,m/10,20,m/10), labels=c("FP","TP","TN","FN"))
hist(pvals24, breaks=20, ylim=c(0,m/3), main="B")
abline(h=(1000-400)/20, col=3)
abline(v = 0.05, col=2)
text(x = c(0.02, 0.02, 0.2, 0.2),y=c(20,100,20,100), labels=c("FP","TP","TN","FN"))

xx <- data.frame(p.vals = pvals24, FDR = p.adjust(pvals24, method = "BH"))

xx <- xx[sample(1:m,floor(m*2/3)),]
xx$FDR2 <- p.adjust(xx$p.vals, method="BH")


par(mfrow=c(1,2))
plot(xx$FDR,xx$FDR2)
abline(c(0,1), col=2)
plot(xx$FDR, xx$FDR-xx$FDR2)
abline(h=0)


length(pvals24)
pvals240 <- c(pvals24, rep(1,300))
hist(pvals240)
xx0 <- data.frame(p.value = pvals240, FDR = p.adjust(pvals240, method = "BH"))

xx0F <- xx0[-(3001:3300),]
nrow(xx0F)
xx0F <- xx0F[sample(1:m,floor(m*2/3)),]
xx0F$FDR2 <- p.adjust(xx0F$p.value, method = "BH")
par(mfrow=c(1,2))
hist(xx0$p.value)
hist(xx0F$p.value)

plot(xx0F$FDR, xx0F$FDR2)
plot(xx0F$FDR2,xx0F$FDR - xx0F$FDR2)
