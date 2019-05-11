# Reference: https://katatrepsis.com/2018/09/24/elisa-analysis-in-r/

#Install the package and load it into the workspace

install.packages("drc")
library(drc)

# For the sake of this example, here is some data from my lab:

dat<-data.frame(Conc=rep(c(0.0001,0.15,0.4,1,2,5),each=3),
                OD=c(2.45,2.51,2.52,2.07,2.15,2,1.59,1.44,1.53,
                     1.05,1.1,1.04,0.661,0.769,0.822,0.56,0.547,0.567))

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=dat)
plot(model1)

# Here are some observed values from new samples
response<-c(2.45,1.67,1.42)
# Estimate the concentration
DOSEx<-ED(model1,response,type="absolute",display=F)
# And here are the estimates including standard errors
DOSEx

# We can add those to our plot
points(y=response,x=DOSEx[,1],col="blue",pch=19,cex=2)

# With error bars
arrows(DOSEx[,1],response,DOSEx[,1]+DOSEx[,2]*1.96,response,length=0.1,angle=90,lwd=3,col="blue")
arrows(DOSEx[,1],response,DOSEx[,1]-DOSEx[,2]*1.96,response,length=0.1,angle=90,lwd=3,col="blue")


## My Code

#first the sample data
#note that field sep might be different based on pre-formatting
cn <- read.csv("Samples2.csv", header=TRUE)
cn[cn < 0] <- 0
#then the standards:
cn_std <- read.csv("STD_FI_Bkgd_vs_ExpCont(2).csv", header=TRUE)

# Create a LL3 model 
fit <- drm(cn_std$Chitinase3_like1_72_FI_Bkgd ~ cn_std$Chitinase3_like1_72_ExpConc ,fct = LL2.3())
summary(fit)
plot(fit)


# Estimate concentration
DOSEx<-ED(fit,cn$Osteopontin_OPN__77_, type="absolute", bound = F)
points(y=cn$Osteopontin_OPN__77_, x=DOSEx[,1],col="blue",pch=19,cex=1)

Results$Osteopontin_OPN__77_ <- DOSEx

