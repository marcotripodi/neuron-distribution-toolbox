set.seed=125

# helper function
read.subsamples <- function(fs) { # keep the same n for each dataset (randomly sub-sample larger sets)

	dat <- list()
    for(i in 1:length(fs)) {
        dat[[i]] <- read.delim(fs[i], header=TRUE)
    }
    n <- unlist(lapply(dat,nrow))
    dat.sub <- NULL
    for(i in 1:length(dat)) {
	    dat.sub <- rbind(dat.sub, cbind(sample=rep(i,min(n)), dat[[i]][sample(n[i],size=min(n)),]))
    }
    dat.sub
}
read.subsamples2 <- function(fs) { # keep all (no sub-sampling to have equal n)
	dat <- NULL
    for(i in 1:length(fs)) {
        dat <- rbind(dat, read.delim(fs[i], header=TRUE))
    }
    dat
}


# read data
setwd("/Users/marco/Documents/RWD/INs coordinates/WT")
# insert the path of your txt files here between the quotation marks separated by a comma
VLfiles <- c(" ")
VL <- read.subsamples2(VLfiles)

BFfiles <- c("")
BF <- read.subsamples2(BFfiles)


GSfiles <- c("")
GS <- read.subsamples2(GSfiles)

TAfiles <- c("")
TA <- read.subsamples2(TAfiles)



GRfiles <- c("")
GR <- read.subsamples2(GRfiles)



#TAKE ONLY POSITIVE VALUES (IPSILATERAL NEURONS)

VL<-VL[VL$x>=0,]
BF<-BF[BF$x>=0,]
GS<-GS[GS$x>=0,]
TA<-TA[TA$x>=0,]
GR<-GR[GR$x>=0,]

#TRANSFORM PIXELS IN MICROMETERS [USING 4X LENS 1PX=2.302 MICRONS]

VL[,'x'] <- 2.302*VL[,'x']
VL[,'y'] <- 2.303*VL[,'y']

BF[,'x'] <- 2.302*BF[,'x']
BF[,'y'] <- 2.302*BF[,'y']

GS[,'x'] <- 2.302*GS[,'x']
GS[,'y'] <- 2.302*GS[,'y']

TA[,'x'] <- 2.302*TA[,'x']
TA[,'y'] <- 2.302*TA[,'y']

GR[,'x'] <- 2.302*GR[,'x']
GR[,'y'] <- 2.302*GR[,'y']



# bin data according to column "t"
bins <- seq(0,150,by=15)
VL <- VL[VL$t>=min(bins) & VL$t <=max(bins),]
BF <- BF[BF$t>=min(bins) & BF$t <=max(bins),]
GS <- GS[GS$t>=min(bins) & GS$t <=max(bins),]
TA <- TA[TA$t>=min(bins) & TA$t <=max(bins),]
GR <- GR[GR$t>=min(bins) & GR$t <=max(bins),]
b.VL <- cut(VL$t, bins, include.lowest=TRUE)
b.BF <- cut(BF$t, bins, include.lowest=TRUE)
b.GS <- cut(GS$t, bins, include.lowest=TRUE)
b.TA <- cut(TA$t, bins, include.lowest=TRUE)
b.GR <- cut(GR$t, bins, include.lowest=TRUE)

quartz(width=10, height=13, pointsize=18)
par(mfrow=c(ceiling((length(bins)-1)/2),2))
for(lev in levels(b.VL)) {
	d.VL <- density(VL[b.VL==lev,"x"])
	d.BF <- density(BF[b.BF==lev,"x"])
	ylims <- range(0,d.VL$y,d.BF$y)
	plot(d.VL, ylim=ylims, xlim=c(0,600), xlab="M-L position", col="red",
	     ylab="INs density", main=sprintf("VL/BF %s",as.character(lev)))
	lines(d.BF, col="green")
	abline(v=0, lty=1, col="black")
	abline(v=median(VL[b.VL==lev,"x"]), lty=3,lwd=3, col="red")
	abline(v=median(BF[b.BF==lev,"x"]), lty=3,lwd=3, col="green")
	legend(x="topright", bty="n", legend=sprintf("P = %.3g", wilcox.test(VL[b.VL==lev,"x"],BF[b.BF==lev,"x"])$p.value))
}


quartz(width=6, height=10, pointsize=10)
par(mfrow=c(ceiling((length(bins)-1)/2),2))
for(lev in levels(b.GS)) {
	d.GS <- density(GS[b.GS==lev,"x"])
	d.TA <- density(TA[b.TA==lev,"x"])
	ylims <- range(0,d.GS$y,d.TA$y)
	plot(d.GS, ylim=ylims, xlim=c(0,600), xlab="M-L position",
	     ylab="INs density", main=sprintf("GS/TA %s",as.character(lev)),col="red")
	lines(d.TA, col="green")
	abline(v=0, lty=1, col="black")
	abline(v=median(GS[b.GS==lev,"x"]), lty=3,lwd=3, col="red")
	abline(v=median(TA[b.TA==lev,"x"]), lty=3,lwd=3, col="green")
	legend(x="topright", bty="n", legend=sprintf("P = %.3g", wilcox.test(GS[b.GS==lev,"x"],TA[b.TA==lev,"x"])$p.value))
}

quartz(width=6, height=10, pointsize=10)
par(mfrow=c(ceiling((length(bins)-1)/2),2))
for(lev in levels(b.GS)) {
	d.GS <- density(GS[b.GS==lev,"x"])
	d.TA <- density(TA[b.TA==lev,"x"])
	d.GR <- density(GR[b.GR==lev,"x"])
	ylims <- range(0,d.GS$y,d.TA$y)
	plot(d.GS, ylim=ylims, xlim=c(0,600), xlab="M-L position",
	     ylab="INs density", main=sprintf("GS/TA %s",as.character(lev)),col="red")
	lines(d.TA, col="green")
	lines(d.GR, col="blue")
	abline(v=0, lty=1, col="black")
	abline(v=median(GS[b.GS==lev,"x"]), lty=3,lwd=3, col="red")
	abline(v=median(TA[b.TA==lev,"x"]), lty=3,lwd=3, col="green")
	abline(v=median(GR[b.GR==lev,"x"]), lty=3,lwd=3, col="blue")
	legend(x="topright", bty="n", legend=sprintf("P = %.3g", wilcox.test(GS[b.GS==lev,"x"],GR[b.GR==lev,"x"])$p.value), GS[b.GS==lev,"x"])
}



quartz(width=10, height=13, pointsize=18)
par(mfrow=c(ceiling((length(bins)-1)/2),2))
for(lev in levels(b.VL)) {
	d.VL <- density(VL[b.VL==lev,"x"])
	d.GS <- density(GS[b.GS==lev,"x"])
	ylims <- range(0,d.VL$y,d.GS$y)
	plot(d.VL, ylim=ylims, xlim=c(0,600), xlab="M-L position",col="green",
	     ylab="INs density", main=sprintf("VL/GS %s",as.character(lev)))
	lines(d.GS, col="red")
	abline(v=0, lty=1, col="black")
	abline(v=median(VL[b.VL==lev,"x"]), lty=3,lwd=3, col="green")
	abline(v=median(GS[b.GS==lev,"x"]), lty=3,lwd=3, col="red")
	legend(x="topright", bty="n", legend=sprintf("P = %.3g", wilcox.test(VL[b.VL==lev,"x"],GS[b.GS==lev,"x"])$p.value))
}

quartz(width=10, height=13, pointsize=18)
par(mfrow=c(ceiling((length(bins)-1)/2),2))
for(lev in levels(b.BF)) {
	d.BF <- density(BF[b.BF==lev,"x"])
	d.TA <- density(TA[b.TA==lev,"x"])
	ylims <- range(0,d.BF$y,d.TA$y)
	plot(d.BF, ylim=ylims, xlim=c(0,600), xlab="M-L position",
	     ylab="INs density", main=sprintf("BF/TA %s",as.character(lev)),col="green")
	lines(d.TA, col="red")
	abline(v=0, lty=1, col="black")
	abline(v=median(BF[b.BF==lev,"x"]), lty=3,lwd=3, col="green")
	abline(v=median(TA[b.TA==lev,"x"]), lty=3,lwd=3, col="red")
	legend(x="topright", bty="n", legend=sprintf("P = %.3g", wilcox.test(GS[b.GS==lev,"x"],TA[b.TA==lev,"x"])$p.value))
}

MEDIAN COMPARISON

# compare difference of medians per bin
quartz(width=10, height=13, pointsize=18)
med.VL <- unlist(lapply(split(VL$x,b.VL),median))
med.BF <- unlist(lapply(split(BF$x,b.BF),median))
barplot(med.BF-med.VL, names.arg=names(med.VL), col="black")

# compare difference of medians per bin
quartz(width=10, height=13, pointsize=18)
med.GS <- unlist(lapply(split(GS$x,b.GS),median))
med.TA <- unlist(lapply(split(TA$x,b.TA),median))
barplot(med.TA-med.GS, names.arg=names(med.GS), col="black")

# compare difference of medians per bin
quartz(width=10, height=13, pointsize=18)
med.GS <- unlist(lapply(split(GS$x,b.GS),median))
med.VL <- unlist(lapply(split(VL$x,b.VL),median))
barplot(med.VL-med.GS, names.arg=names(med.GS), col="black")

# compare difference of medians per bin
quartz(width=10, height=13, pointsize=18)
med.TA <- unlist(lapply(split(TA$x,b.TA),median))
med.BF <- unlist(lapply(split(BF$x,b.BF),median))
barplot(med.TA-med.BF, names.arg=names(med.TA), col="black")

MAX COMPARISON

barplot(c(max(TA$y)-max(GS$y),max(BF$y)-max(VL$y)))

