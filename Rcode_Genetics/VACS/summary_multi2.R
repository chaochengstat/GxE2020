setwd("/gpfs/ycga/project/wang_zuoheng/cc2747/VACS/spline200117")
library("qqman")
library(doMC)
library(doRNG)
library("snpStats")
library("car")
library("caret")
library(sandwich)
library(lmtest)
library("GEint")
library("foreach")


registerDoMC(20)
registerDoRNG(1)

files=dir()
split=strsplit(files, split=" ")
loc=rep(0,length(files))
for (i in (1:length(loc))) {
  if (split[[i]][1]=="The" & split[[i]][3]=="manhattan") {
    loc[i]=1
  }
}
files=files[which(loc==1)]#[seq(from=1,to=340,by=10)]#[c(1,10,100,200,300,342)]

load(files[1]) 
manhattan.data=manhattan
for (i in (2: length(files))) {
  load(files[i]) 
  x=manhattan
  manhattan.data=rbind(manhattan.data,x)
}


#manhattan.data=foreach(i=1:length(files),.combine='rbind') %dopar% {
#  print(i)
#  load(files[i]) 
#  x
#}

print(dim(manhattan.data))
##### parameters for the plot
width=1300;height=800;cex=2

print("smoking")

colnames(manhattan.data)[-c(1:3)]=paste("p",1:8,sep="")

##### original data
### (a) reverse test

labels_a=c("(a) Reverse test for smoking based on all samples (N=2161 with 10,592,437 SNPs)",
           "(a) Reverse test for alcohol based on all samples (N=2161 with 10,592,437 SNPs)")
labels_c=c("(c) Logistic test for smoking based on all samples (N=2161 with 10,592,437 SNPs)",
           "(c) Logistic test for alcohol based on all samples (N=2161 with 10,592,437 SNPs)")


output = foreach (i=(1:2),.combine='rbind') %dopar% {
  
  j=2*i;m=2*i-1
  mname1=paste("p",j,sep="");mname2=paste("p",m,sep="")
  min.p=min(min(manhattan.data[,mname1],na.rm=T),min(manhattan.data[,mname2],na.rm=T))[1]
  min.p=min(min.p,4*10^(-8))
  ylim=c(0,-log10(min.p)+0.25)
  manhattan.data1=na.omit(manhattan.data[,c("chromosome","position","snp.name",mname1)])
  manhattan.data2=na.omit(manhattan.data[,c("chromosome","position","snp.name",mname2)])
  
  
  chisq1 = qchisq(manhattan.data1[,mname1],1,lower.tail=FALSE);
  chisq2 = qchisq(manhattan.data2[,mname2],1,lower.tail=FALSE);
  options(bitmapType='cairo')
  lam1=round(median(chisq1) / qchisq(0.5,1),3)
  lam2=round(median(chisq2) / qchisq(0.5,1),3)

  if ((lam1>1.05) |(lam1 < 0.95)) {
    manhattan.data1[,mname1] = pchisq(chisq1/lam1, df = 1, lower = F)
  }
  
  if ((lam2>1.05) |(lam2 < 0.95)) {
    manhattan.data2[,mname2] = pchisq(chisq2/lam2, df = 1, lower = F)
  }

  png(paste("double_",mname1,".png",sep=""),width=width,height=height)
  
  par(mar=c(6,6,6,6))
  layout(matrix(c(1,1,1,1,1,2,2,3,3,3,3,3,4,4),ncol=7,byrow=T))
  qqman::manhattan(x=manhattan.data1,chr="chromosome",bp="position",p=mname1,ylim=ylim,
                   snp="snp.name",cex.axis=cex,cex.lab=cex,cex.main=cex,main=labels_a[i])
  qqman::qq(manhattan.data1[,mname1],cex.axis=cex,cex.lab=cex,cex.main=cex,main="(b) QQ plot for reverse test")
  #t=bquote(bold(lambda == .(lam1)))
  #y=-log10(min(as.numeric(manhattan.data1[,mname1])))*0.8
  #text(1,y,t,cex=2)
  
  qqman::manhattan(x=manhattan.data2,chr="chromosome",bp="position",p=mname2,ylim=ylim,
                   snp="snp.name",cex.axis=cex,cex.lab=cex,cex.main=cex,main=labels_c[i])
  qqman::qq(manhattan.data2[,mname2],cex.axis=cex,cex.lab=cex,cex.main=cex,main="(d) QQ plot for logistic test")
  #t=bquote(bold(lambda == .(lam2)))
  #y=-log10(min(as.numeric(manhattan.data2[,mname2])))*0.8
  #text(1,y,t,cex=2)
  
  
  dev.off()
  1
}

loc=which(apply(manhattan.data,1,function(x) min(as.numeric(x[-c(1:4)])))<1e-6)
selectedsnp=manhattan.data[loc,]
write.csv(selectedsnp,"selectedsnp.csv")