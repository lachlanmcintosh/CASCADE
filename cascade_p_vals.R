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
data$cn2 <- data$cn.float * (1-data$maf.float)
data$var1 <- sapply(data$segment,function(x) sd(data[which(data$segment == x), "cn1"])^2)
data$var2 <- sapply(data$segment,function(x) sd(data[which(data$segment == x), "cn2"])^2)
max(data$var1,na.rm=TRUE)
max(data$var2,na.rm=TRUE)

min(data$var1,na.rm=TRUE)
min(data$var2,na.rm=TRUE)

mean(data$var1,na.rm=TRUE)
mean(data$var2,na.rm=TRUE)

unique = which(!duplicated(data$segment))
pooled_var <- sum((data[unique,"n"] -rep(1,nrow(data[unique,])))*data[unique,"var1"],na.rm=TRUE)/(nrow(data) - max(data$segment))

data$cnA <- apply(data[,c("cn1","cn2")],1,min)
data$cnB <- apply(data[,c("cn1","cn2")],1,max)

data$cn_roundedA <- round(data$cnA,0)
data$cn_roundedB <- round(data$cnB,0)

data$tA <- (data$cn_roundedA - data$cnA)/(pooled_var)^0.5
data$tB <- (data$cn_roundedB - data$cnB)/(pooled_var)^0.5
data$pA <- pt(data$tA,rep(1,nrow(data)))
data$pB <- pt(data$tB,rep(1,nrow(data)))
data$p <- apply(data[,c("pA","pB")],1,min)


#ggplot(data)+geom_point(aes(x=loc,y=cnA),col=data$pA)+geom_point(aes(x=loc,y=cnB),col=data$pB)+facet_wrap(~sample+chr,ncol=4)
un <- unique(data$sample)

data$chr <- as.factor(data$chr)
ggplot(data[which(data$sample == un[1]),])+
    geom_point(aes(x=loc,y=cnB,col=pB))+
    geom_point(aes(x=loc,y=cnA, col=pA))+
    facet_wrap(~chr,scales="free",ncol=3)

max(data$pA)
min(data$pA)
max(data$pB)
min(data$pB)


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
 
summ6 <- data[unique(data$segment),c("segment","chr","sample","var1","var2")]

summ <- merge(x=merge(x=merge(x=merge(x=merge(x=summ,y=summ2,by=1),y=summ3,by=1),y=summ4,by=1),y=summ5,by=1),y=summ6,by=1)


summ$cnA_rounded <- round(summ$cnA)
summ$cnB_rounded <- round(summ$cnB)
summ$tA <- (summ$cnA_rounded - summ$cnA)/(summ$n/summ$var1)^0.5
summ$tB <- (summ$cnB_rounded - summ$cnB)/(summ$n/summ$var2)^0.5

summ$pA <- pt(summ$tA,summ$n)
summ$pB <- pt(summ$tB,summ$n)

summ$chr <- as.factor(summ$chr)
ggplot(summ[which(summ$sample == un[1]),])+
  geom_segment(aes(x=minx,y=cnA,xend=maxx,yend=cnA,col=pA))+
  geom_segment(aes(x=minx,y=cnB,xend=maxx,yend=cnB,col=pB))+
  facet_wrap(~chr,scales="free",ncol=3)
