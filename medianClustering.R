# read data
setwd("/Users/marco/Documents/RWD/INs coordinates/WT")

# global parameters

# use these for correleation segment by segment
# bins <- seq(0,150,by=15)
# binnms <- 1:10

# use the following for correleation all segments projected on one plane
bins <- c(0,150)
binnms <- 1


xmin <- 0
tmax <- max(bins)


# Insert files here

# Insert name of the distinct distribution that you intend to analyse
types <- c("VL","GS","BF","TA","GR")

# Insert the text files containing the coordinates of your neuronal populations between quotes. One file per reconstruction, as many files as you want. Here you'll merge all files in one large merged file
VLfiles <- c("file1.text, file2.txt ")
BFfiles <- c(" file1.text, file2.txt ")
GSfiles <- c(" file1.text, file2.txt ")
TAfiles <- c("file1.text, file2.txt ")
GRfiles <- c(" file1.text, file2.txt ")

# load data and calculate per mouse median for each bin
dat <- NULL
for(type in types) {
    for(infile in get(paste(type,"files",sep=""))) {
	    tmp <- read.delim(infile, header=TRUE)
		tmp <- tmp[tmp$x>=xmin & tmp$t<=tmax, ]
		tmp2 <- tapply(tmp$x, findInterval(tmp$t, bins, rightmost.closed=TRUE), median)
		i <- binnms[as.integer(names(tmp2))]
		n <- length(i)
		dat <- rbind(dat, data.frame(file=rep(infile,n), muscle=rep(type,n), zbin=i, med=tmp2))
	}
}

tab <- list()
pvals <- list()
for(b in binnms) {
    tab[[b]] <- matrix(NA, nrow=length(types), ncol=length(types), dimnames=list(types, types))
	pvals[[b]] <- matrix(NA, nrow=length(types), ncol=length(types), dimnames=list(types, types))
    for(t1 in types) {
    	for(t2 in types) {
	    	tab[[b]][t1,t2] <- abs(mean(dat[dat$muscle==t1 & dat$zbin==b,'med']) - mean(dat[dat$muscle==t2 & dat$zbin==b,'med']))
			if(length(dat[dat$muscle==t1 & dat$zbin==b,'med'])>1 && length(dat[dat$muscle==t2 & dat$zbin==b,'med'])>1)
    	    	pvals[[b]][t1,t2] <- t.test(dat[dat$muscle==t1 & dat$zbin==b,'med'], dat[dat$muscle==t2 & dat$zbin==b,'med'])$p.value
		}
	}
}

library(RColorBrewer)
plotcols <- rev(brewer.pal(9, "Greys"))

for(b in binnms) {
	quartz()
    image(t(tab[[b]]), col=plotcols, axes=FALSE, main=b)
    xs <- seq(0,1,length=length(types))
    axis(1,at=xs,label=types)
    axis(2,at=xs,label=types)
    box()
    for(i in 1:length(types)) {
    	for(j in 1:length(types)) {
    		if(i!=j) {
    		   text(x=xs[i], y=xs[j], adj=c(.5,.5),
    			    label=sprintf("P=%.3g",pvals[[b]][i,j]), col="red")
    	    }
    	}
    }
}















