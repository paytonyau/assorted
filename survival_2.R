## http://www.biotrainee.com/thread-929-1-1.html

A=data.frame()

#�������
x=c(2, -2.5, 3.5, 4, 4, -5, 6, -6, 7, -7, 8, -9, 10.5, 12.5, 19, 2.5, 5, 
    7, -8.5, 9, -10, 11, -11, 12, 13, -14, 15, -16, 17, -18, 19, -20, 21, 24, 32)

#��survtС��0(��������ɾʧ)�������censor��ֵΪ1�����򣬸�ֵΪ0
for (i in 1:length(x))
{
  A[i,1]=x[i]
  if(A[i,1]<0)
    A[i,2]=1
  else
    A[i,2]=0
}

#��Ϊhigh-wbc��low-wbc����
A[,3] <- c(rep("high",15),rep("low",20))
A[,1] <- abs(x)

#��������
names(A) <- c("survt","censor","wbc")


# �����������
library(survival)
Surv(A$survt,A$censor==1)

# ����KM��������(Kaplan�CMeier estimator�˻����޷�)
y <- Surv(A$survt,A$censor==1)
kmfit1 <- survfit(y~1)
summary(kmfit1)
plot(kmfit1)

# ����wbc�������KM��������
kmfit2 <- survfit(y~A$wbc)
plot(kmfit2, lty = c('solid', 'dashed'), col=c('black','blue'),
     xlab='survival time in days',ylab='survival probabilities')
legend('topright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))
summary(kmfit2)

# ����������,rho=0Ϊlog-rank����Mantel Haenszel����rho=1ΪWilcoxon��
survdiff(Surv(survt,censor)~wbc, data=A,rho = 0)
survdiff(Surv(survt,censor)~wbc, data=A,rho = 1)

#���ò�ͬ�����Ĳ�ͬ������ѡ��ͬ�ļ��鷽��
survreg(Surv(survt,censor)~wbc, data=A,dist="weibull")
survreg(Surv(survt,censor)~wbc, data=A,dist="logistic")
survreg(Surv(survt,censor)~wbc, data=A,dist="lognormal")

#-logS(p)������ʱ��t��ɢ��ͼ
kmfit3=kmfit2
kmfit3$surv=-log(kmfit3$surv)
b1=data.frame(kmfit3$time[1:12],kmfit3$surv[1:12])
b2=data.frame(kmfit3$time[13:31],kmfit3$surv[13:31])
plot(b1,type="b",col="black",main="�������溯���ĸ�������")
lines(b2,type="b",col="blue")
legend('bottomright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))

#log��-logS(t)����log(t)��ɢ��ͼ
kmfit4=kmfit3
kmfit4$surv=log(kmfit4$surv)
C1=data.frame(kmfit4$time[1:12],kmfit4$surv[1:12])
C2=data.frame(kmfit4$time[13:31],kmfit4$surv[13:31])
plot(C1,type="b",col="black",main="�������溯���ĸ��������Ķ���")
lines(C2,type="b",col="blue")
legend('bottomright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))


# ��word�ļ���Ĺ�ʽ--��n.risk��ֵ��/��n.risk*time��ֵ����
# ����wbc�ֱ�Ϊhigh��low���������պ���ͼ
summary(kmfit2)
B =data.frame(kmfit2$time[1:12],kmfit2$n.risk[1:12])
B$rr= 0
for (i in 2:12)
{
  B[i,3]=(B[i-1,2]-B[i,2])/(B[i-1,2]*(B[i,1]-B[i-1,1]))
}

plot(B$rr~B$kmfit2.time.1.12.,xlab = "t",ylab = "������",col="black",type="l",main="�������պ���")
legend('topright', c('wbc high','wbc low'), lty=c('solid','dashed'),
       col=c('black','blue'))

C =data.frame(kmfit2$time[13:31],kmfit2$n.risk[13:31])
C$rr= 0
for (i in 2:19)
{
  C[i,3]=(C[i-1,2]-C[i,2])/(C[i-1,2]*(C[i,1]-C[i-1,1]))
}

lines(C$rr~C$kmfit2.time.13.31.,col="blue",type="l")