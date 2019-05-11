## http://www.biotrainee.com/thread-929-1-1.html

A=data.frame()

#输入变量
x=c(2, -2.5, 3.5, 4, 4, -5, 6, -6, 7, -7, 8, -9, 10.5, 12.5, 19, 2.5, 5, 
    7, -8.5, 9, -10, 11, -11, 12, 13, -14, 15, -16, 17, -18, 19, -20, 21, 24, 32)

#若survt小于0(代表数据删失)，则变量censor赋值为1，否则，赋值为0
for (i in 1:length(x))
{
  A[i,1]=x[i]
  if(A[i,1]<0)
    A[i,2]=1
  else
    A[i,2]=0
}

#分为high-wbc和low-wbc两组
A[,3] <- c(rep("high",15),rep("low",20))
A[,1] <- abs(x)

#变量命名
names(A) <- c("survt","censor","wbc")


# 建立生存对象
library(survival)
Surv(A$survt,A$censor==1)

# 估计KM生存曲线(Kaplan–Meier estimator乘积极限法)
y <- Surv(A$survt,A$censor==1)
kmfit1 <- survfit(y~1)
summary(kmfit1)
plot(kmfit1)

# 根据wbc分组估计KM生存曲线
kmfit2 <- survfit(y~A$wbc)
plot(kmfit2, lty = c('solid', 'dashed'), col=c('black','blue'),
     xlab='survival time in days',ylab='survival probabilities')
legend('topright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))
summary(kmfit2)

# 检验显著性,rho=0为log-rank法或Mantel Haenszel法，rho=1为Wilcoxon法
survdiff(Surv(survt,censor)~wbc, data=A,rho = 0)
survdiff(Surv(survt,censor)~wbc, data=A,rho = 1)

#设置不同函数的不同参数，选择不同的检验方法
survreg(Surv(survt,censor)~wbc, data=A,dist="weibull")
survreg(Surv(survt,censor)~wbc, data=A,dist="logistic")
survreg(Surv(survt,censor)~wbc, data=A,dist="lognormal")

#-logS(p)对生存时间t的散点图
kmfit3=kmfit2
kmfit3$surv=-log(kmfit3$surv)
b1=data.frame(kmfit3$time[1:12],kmfit3$surv[1:12])
b2=data.frame(kmfit3$time[13:31],kmfit3$surv[13:31])
plot(b1,type="b",col="black",main="估计生存函数的负数对数")
lines(b2,type="b",col="blue")
legend('bottomright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))

#log（-logS(t)）对log(t)的散点图
kmfit4=kmfit3
kmfit4$surv=log(kmfit4$surv)
C1=data.frame(kmfit4$time[1:12],kmfit4$surv[1:12])
C2=data.frame(kmfit4$time[13:31],kmfit4$surv[13:31])
plot(C1,type="b",col="black",main="估计生存函数的负数对数的对数")
lines(C2,type="b",col="blue")
legend('bottomright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))


# 由word文件里的公式--（n.risk差值）/（n.risk*time差值），
# 画出wbc分别为high和low的死亡风险函数图
summary(kmfit2)
B =data.frame(kmfit2$time[1:12],kmfit2$n.risk[1:12])
B$rr= 0
for (i in 2:12)
{
  B[i,3]=(B[i-1,2]-B[i,2])/(B[i-1,2]*(B[i,1]-B[i-1,1]))
}

plot(B$rr~B$kmfit2.time.1.12.,xlab = "t",ylab = "死亡率",col="black",type="l",main="死亡风险函数")
legend('topright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))

C =data.frame(kmfit2$time[13:31],kmfit2$n.risk[13:31])
C$rr= 0
for (i in 2:19)
{
  C[i,3]=(C[i-1,2]-C[i,2])/(C[i-1,2]*(C[i,1]-C[i-1,1]))
}

lines(C$rr~C$kmfit2.time.13.31.,col="blue",type="l")
