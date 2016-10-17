rm(list=ls())
data <- read.csv("/Users/lmcintosh/Documents/CASCADE/CA0053_tumour_est_subclonal.txt",sep="\t")

data$cn.float <- as.numeric(data$cn.float)
data$maf.float <- as.numeric(data$maf.float)
#data$cn.float <- as.numeric(data$cn.float)
#data$cn.float <- as.numeric(data$cn.float)
data$chr <- sapply(data$mutation, function(x) gsub( ";.*$", "", x ))
data$loc <- sapply(as.character(data$mutation), function(x) as.numeric(gsub("^.*;([0-9]+).*$", "\\1", x)))

data$segment <- rep(0,nrow(data))
current <- 1
for(i in 1:nrow(data)){
  if(i>1 && (data[i,"cn.float"] != data[i-1,"cn.float"] || data[i,"chr"] != data[i-1,"chr"] || data[i,"sample"] != data[i-1,"sample"])){
    current = current + 1
  }
  data[i,"segment"] <- current 
}
length(unique(data$segment))
length(unique(data$chr))+length(unique(data$cn.float))
unique(data$sample)

data$n <- sapply(data$segment,function(x) sum(data$segment == x))
data$cn1 <- data$cn.float * data$maf.float
#data$cn2 <- data$cn.float * (1-data$maf.float)

start = 0.5
end = 2
by = 0.0001
range = seq(start,end,by)
#best_pooled_var = infty()
pooled_var_vec <- c()
for(i in range){
  #data2 <- c(((data[, "cn1"]+i/2)%%i),((data[, "cn2"]+i/2)%%i))/i -0.5
  data2 <- ((data[, "cn1"]+i/2)%%i)/i -0.5
  
  #ggplot()+geom_histogram(aes(data2))
  pooled_var <- (sum((data2)^2)/(length(data2)-1))
  print(i)
  if(i==range[1] || pooled_var < best_pooled_var){
    best = i
    best_pooled_var = pooled_var
  }
  pooled_var_vec <- c(pooled_var_vec,pooled_var)
}
i = best
data2<- (data[, "cn1"]+i/2)%%i
pooled_var <- sd(data2)^2
library(ggplot2)
ggplot()+geom_histogram(aes(data2))

pooled_var_norm <- (sum((data2-0.5)^2)/(length(data2)-1))


data$var1 <- sapply(data$segment,function(x) sd(data[which(data$segment == x), "cn1"])^2)
unique = which(!duplicated(data$segment))
pooled_var1 <- 0.5*sum((data[unique,"n"] -rep(1,nrow(data[unique,])))*data[unique,"var1"],na.rm=TRUE)/(nrow(data) - max(data$segment))
pooled_var_crap <- pooled_var1
changes <- pooled_var_vec[1:(length(pooled_var_vec)-1)] - pooled_var_vec[2:(length(pooled_var_vec))]
changes2 <- changes[1:(length(changes)-1)] - changes[2:(length(changes))]

library(ggplot2)
ggplot()+geom_point(aes(x=range,y=pooled_var_vec)) + geom_hline(yintercept = pooled_var_crap)+ geom_hline(yintercept = pooled_var)+ geom_hline(yintercept = pooled_var_norm)#+ylim(0,1)
#ggplot()+geom_point(aes(x=range,y=pooled_var_vec)) + geom_hline(yintercept = pooled_var)+ylim(0,0.01)

pooled_var <- pooled_var_norm

ggplot()+geom_point(aes(x=range[1:(length(range)-1)],y=changes))
ggplot()+stat_smooth(aes(x=50:148/100,y=changes[1:99],span=0.05))
ggplot()+geom_point(aes(x=50:147/100,y=abs(changes2[1:98])))
ggplot()+stat_smooth(aes(x=50:147/100,y=changes2[1:98],span=0.05))


# data$var1 <- sapply(data$segment,function(x) sd(data[which(data$segment == x), "cn1"])^2)
# 
# 
# unique = which(!duplicated(data$segment))
# pooled_var <- sum((data[unique,"n"] -rep(1,nrow(data[unique,])))*data[unique,"var1"],na.rm=TRUE)/(nrow(data) - max(data$segment))


data$cn_rounded <- round(data$cn1,0)

data$t <- (data$cn_rounded - data$cn1)/(pooled_var)^0.5
data$p <- 2*pt(-abs(data$t),rep(1,nrow(data)))


#ggplot(data)+geom_point(aes(x=loc,y=cnA),col=data$pA)+geom_point(aes(x=loc,y=cnB),col=data$pB)+facet_wrap(~sample+chr,ncol=4)
un <- unique(data$sample)

data$chr <- as.factor(data$chr)
ggplot(data[which(data$sample == un[1]),])+
  geom_point(aes(x=loc,y=cn1,col=p))+
  facet_wrap(~chr,scales="free",ncol=3)

max(data$p)
min(data$p)


# we need to aggregate data here:
attach(data)
summ <- aggregate(data[,c("segment","loc")],by=list(segment),FUN=min)
summ <- summ[,c("segment","loc")]
colnames(summ) <- c("segment","minx")

summ2 <- aggregate(data[,c("segment","loc")],by=list(segment),FUN=max)
summ2 <- summ2[,c("segment","loc")]
colnames(summ2) <- c("segment","maxx")

summ3 <- aggregate(data[,c("segment","cnA")],by=list(segment),FUN=mean)
summ3 <- summ3[,c("segment","cnA")]

summ4 <- aggregate(data[,c("segment","cnB")],by=list(segment),FUN=mean)
summ4 <- summ4[,c("segment","cnB")]

summ5 <- aggregate(data[,c("segment","n")],by=list(segment),FUN=max)
summ5 <- summ5[,c("segment","n")]

summ6 <- data[unique(data$segment),c("segment","chr","sample")]

summ <- merge(x=merge(x=merge(x=merge(x=merge(x=summ,y=summ2,by=1),y=summ3,by=1),y=summ4,by=1),y=summ5,by=1),y=summ6,by=1)


summ$cnA_rounded <- round(summ$cnA)
summ$cnB_rounded <- round(summ$cnB)
summ$tA <- (summ$cnA_rounded - summ$cnA)/(summ$n/pooled_var)^0.5
summ$tB <- (summ$cnB_rounded - summ$cnB)/(summ$n/pooled_var)^0.5

summ$pA <- 2*pt(-abs(summ$tA),summ$n)
summ$pB <- 2*pt(-abs(summ$tB),summ$n)

summ$chr <- as.factor(summ$chr)
ggplot(summ[which(summ$sample == un[1]),])+
  geom_segment(aes(x=minx,y=cnA,xend=maxx,yend=cnA,col=pA))+
  geom_segment(aes(x=minx,y=cnB,xend=maxx,yend=cnB,col=pB))+
  facet_wrap(~chr,scales="free",ncol=3)

facsum(data$p<0.1,na.rm=TRUE)/nrow(data)

# ok now for a visualisation
ggplot()+geom_point()

#should not only do p-values for individual points, should also do p-values for segments.

