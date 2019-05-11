## Adjust P-values for Multiple Comparisons
## Surce http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

require(graphics)

df1 <- read.csv("C:/PYTHON/plateB.csv", header=TRUE)

round(df1, 3)

## Reason
##https://stackoverflow.com/questions/12384071/how-to-coerce-a-list-object-to-type-double
B <- as.numeric(unlist(df1))

round(p.adjust(B), 3) ## 3 sig.fig, chack the data only
round(p.adjust(B, "BH"), 3)

## or all of them at once (dropping the "fdr" alias):
p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
p.adj    <- sapply(p.adjust.M, function(meth) p.adjust(B, meth))
p.adj.60 <- sapply(p.adjust.M, function(meth) p.adjust(B, meth, n = 60))
stopifnot(identical(p.adj[,"none"], B), p.adj <= p.adj.60)
round(p.adj, 3)
## or a bit nicer:

C <- noquote(apply(p.adj, 2, format.pval, digits = 3))

write.csv(C, "C:/PYTHON/pE.csv")

## and a graphic:
matplot(p, p.adj, ylab="p.adjust(p, meth)", type = "l", asp = 1, lty = 1:6,
        main = "P-value adjustments")
legend(0.7, 0.6, p.adjust.M, col = 1:6, lty = 1:6)

## Can work with NA's:
pN <- p; iN <- c(46, 47); pN[iN] <- NA
pN.a <- sapply(p.adjust.M, function(meth) p.adjust(pN, meth))
## The smallest 20 P-values all affected by the NA's :
round((pN.a / p.adj)[1:20, ] , 4)
