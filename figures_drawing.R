library("scales")
library("ggplot2")
library("vioplot")
library("RColorBrewer") 
library("wesanderson")
library("plyr")
####Homework Lab Project2 DDL2018-10-02
##################################################
## Exercise 01 ###################################
##################################################
######read data from file
moviedata = read.csv("moviedata.csv", header=TRUE, sep=",", dec = ".", na.strings ="NA")

##################Deal with NA value
moviedata=na.omit(moviedata)
###create factor variable in the moviedata dataframe
moviedata$sequelcat <- factor(moviedata$dummy_sequel, levels = c(0, 1),
                              labels = c("ORIGINAL", "SEQUEL"))
moviedata$year <- factor(moviedata$year, levels = c(2014, 2015),labels = c("Y2014", "Y2015"))
##################################################
## Exercise 02 ###################################
##################################################
############Draw a boxplot of ratings and year
boxplot(ratings~year, data=moviedata, main="Movie Ratings in 2014,2015", xlab="Year", ylab="Rate Point")
############Make a rate column and violin plot
rate_low <- moviedata$gross[moviedata$ratings<6]
rate_high <- moviedata$gross[moviedata$ratings>8]
rate_mid <- moviedata$gross[moviedata$ratings>=6&moviedata$ratings<=8]
vioplot(rate_low, rate_mid , rate_high, names=c("Rate<6", "Rate6-8", "Rate>8"),col="gold")
title("Violin Plots of Rate and Gross",xlab="Rate",ylab="Gross")
#########Draw barplot of count of rate
counts <- table(as.integer(moviedata$ratings))
barplot(counts, main="movie rate", xlab="movie rate")
##########Draw a scatter plot of rate and count
p <- qplot(gross, budget, data=moviedata)
p + xlab('budget') + ylab('gross')+scale_y_continuous(labels = dollar)+scale_x_continuous(labels = dollar)
winsor <- function(x, multiplier) {
  if(length(multiplier) != 1 || multiplier <= 0) {
    stop("bad value for 'multiplier'")????}
  
  quartile1 = summary(x)[2] # Calculate lower quartile
  quartile3 = summary(x)[5] # Calculate upper quartile
  iqrange = IQR(x) # Calculate interquartile range
  
  y <- x
  boundary1 = quartile1 - (iqrange * multiplier)
  boundary2 = quartile3 + (iqrange * multiplier)
  
  y[ y < boundary1 ] <- boundary1
  y[ y > boundary2 ] <- boundary2
  
  y
}
moviesubset=subset(moviedata)
likes_winsoring=winsor(moviesubset$likes,1.5)
##################################################
## Exercise 03 ###################################
##################################################
####t test
t_test=t.test(moviedata$gross,moviedata$budget)
####anova
A=factor(c(rep(1:6,each=31),1))
anova=aov(gross~A,moviedata)
#####HSD test
hsd=TukeyHSD(anova)
###########wilcox_test
wilcox.test(moviedata$gross, moviedata$budget, alternative="less",
            exact=FALSE,correct=FALSE, conf.int=TRUE)
##################################################
## Exercise 04 ###################################
##################################################
#####create new variable
moviedata$RateStatus[moviedata$ratings >= 0 & moviedata$ratings <= 6] <- "Rate low"
moviedata$RateStatus[moviedata$ratings >= 6 & moviedata$ratings <= 10] <- "Rate high"
moviedata$RateStatus <- as.factor(moviedata$RateStatus)
##################################################
## Exercise 05 ###################################
##################################################
###########aggregate data
ddply(moviedata, c("year"), summarise, N = length(movie), gross_avg=mean(gross), budget_avg=mean(budget), rate_ave=mean(ratings))
