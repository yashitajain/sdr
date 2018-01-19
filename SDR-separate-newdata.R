#SDR For each data type separately 
set.seed(1)
#load("screeningdata400.rda")
load('new_data_anova.rda')
Z=SDRvariable
n=58
p=c(1000,1000,422,162)
d=3
s=4
y=8
#standardization
library(corpcor)

st=list()
for (i in 1:s){
        st[[i]]=cov.shrink(Z[[i]])
}

for (i in 1:s){
        st[[i]]=mpower(st[[i]],-0.5)
}

for (i in 1:s){
        Z[[i]]=(Z[[i]]-mean(Z[[i]]))%*%st[[i]]
}

outcome<- factor(c(rep("1", 5), rep("2", 6),rep("3",7),rep("4",6),rep("5",10),rep("6",9),rep("7",7),rep("8",8)))
outcome=data.frame(outcome)
for (i in 1:s){
        
        Z[[i]]=data.frame(outcome,Z[[i]]) 
        Z[[i]]$outcome=as.factor(Z[[i]]$outcome)
}

cp=Z[[1]]
#mirna
m=matrix(0,p[1],p[1])
for (i in 1:y){
        q=(cp[which (cp$outcome==i),])
        avg=(aggregate(q[,2:1001],list(q$outcome),mean))
        m=m+(t(avg[,-1])%*%as.matrix(avg[,-1]))*sum(with(cp,outcome==i))  
}
m=m/n
decom=eigen(m)
gamma1=decom$vectors[,1:d]


#mrna
mrna=Z[[2]]
m=matrix(0,p[2],p[2])
for (i in 1:y){
        q=(mrna[which (mrna$outcome==i),])
        avg=(aggregate(q[,2:1001],list(q$outcome),mean))
        m=m+(t(avg[,-1])%*%as.matrix(avg[,-1]))*sum(with(mrna,outcome==i))  
}
m=m/n
decom=eigen(m)
gamma2=decom$vectors[,1:d]

mirna=Z[[3]]
m=matrix(0,p[3],p[3])
for (i in 1:y){
        q=(mirna[which (mirna$outcome==i),])
        avg=(aggregate(q[,2:423],list(q$outcome),mean))
        m=m+(t(avg[,-1])%*%as.matrix(avg[,-1]))*sum(with(mirna,outcome==i))  
}
m=m/n
decom=eigen(m)
gamma3=decom$vectors[,1:d]


mirna=Z[[4]]
m=matrix(0,p[4],p[4])
for (i in 1:y){
        q=(mirna[which (mirna$outcome==i),])
        avg=(aggregate(q[,2:163],list(q$outcome),mean))
        m=m+(t(avg[,-1])%*%as.matrix(avg[,-1]))*sum(with(mirna,outcome==i))  
}
m=m/n
decom=eigen(m)
gamma4=decom$vectors[,1:d]

SDR1=as.matrix(Z[[1]][,-1])%*%gamma1
SDR2=as.matrix(Z[[2]][,-1])%*%gamma2
SDR3=as.matrix(Z[[3]][,-1])%*%gamma3
SDR4=as.matrix(Z[[4]][,-1])%*%gamma4




#randomforest
#mirna
sdr=cbind(SDR1,SDR2,SDR3,SDR4,outcome)
colnames(sdr)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','outcome')

#V13','V14','V15','V16'
SDR1=data.frame(SDR1,outcome) 
SDR1$outcome=as.factor(SDR1$outcome)


SDR2=data.frame(SDR2,outcome) 
SDR2$outcome=as.factor(SDR2$outcome)

SDR3=data.frame(SDR3,outcome) 
SDR3$outcome=as.factor(SDR3$outcome)

SDR4=data.frame(SDR4,outcome) 
SDR4$outcome=as.factor(SDR4$outcome)

library(randomForest)
errorset1=c()
rf<-function(x){
pred=NULL
for (i in 1:n){
        test=x[i,]
        train=x[-i,]
        modelfit=randomForest(train$outcome~.,train)
        prediction=predict(modelfit,test)
        pred=rbind(pred,prediction)}
return (pred)
        }




pred0=rf(sdr)
pred1=rf(SDR1)
pred2=rf(SDR2)
pred3=rf(SDR3)
pred4=rf(SDR4)

c<- factor(c(rep("1", 5), rep("2", 6),rep("3",7),rep("4",6),rep("5",10),rep("6",9),rep("7",7),rep("8",8)))

cm=table(pred0,c)
right=diag(cm)
n=sum(cm)
acc=sum(right)/n
acc*100


cm=table(pred1,c)
right=diag(cm)
n=sum(cm)
acc=sum(right)/n
acc*100


cm=table(pred2,c)
right=diag(cm)
n=sum(cm)
acc=sum(right)/n
acc*100


cm=table(pred3,c)
right=diag(cm)
n=sum(cm)
acc=sum(right)/n
acc*100


cm=table(pred4,c)
right=diag(cm)
n=sum(cm)
acc=sum(right)/n
acc*100


