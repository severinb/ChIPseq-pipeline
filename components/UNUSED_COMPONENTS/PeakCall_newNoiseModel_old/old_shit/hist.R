library(Hmisc)

hx <- rnorm(10000000,m=0,sd=1)
d = read.csv(file="outzvals",sep="\t")
x=hist(d[,7],xlim=c(-10,10),col="red",xlab="z-value",breaks=500);

jpeg("hist_zvals.jpg")
  plot(x$mids,x$density,type="s",log="y",xlim=c(-10,10),xlab="Z-values",ylab="Density (log)",main="Histogram: z-values")
  lines(density(hx), type="l", lwd=2, lty=1,col="red")
  n=5
  minor.tick(nx=n, tick.ratio=1)
  #lines( c(2.5,2.5),c(0.0000000000001,1) )
  

dev.off()


#system('R CMD BATCH hist.R');

