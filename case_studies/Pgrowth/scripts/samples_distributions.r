# samples as distributions (and vice-versa)
rm(list=ls())

# parametric distributions (pdfs)
# e.g. normal, uniform

# functional form

x <- seq(-4,4,0.1)
pdfx <- dnorm(x,mean=0,sd=1)

par(mfrow=c(1,1))
plot(x,pdfx,type="b", pch=19)

##############################################################

sampx <- rnorm(500,mean=0,sd=1)
  
# plots  

plot(sampx,type="p", pch=19)

##############################################################

hist(sampx, breaks=10)

##############################################################

ksd <- density(sampx, bw=1, kernel="triangular")
plot(ksd)

