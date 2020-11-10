###########################################
# Functions
# srun --pty --x11 --mem=10000 -p interactive bash
###########################################

trans.yj=function(x) {
  s = x + 1
  m = summary(p1 <- powerTransform(s ~ 1))
  bcPower(s,lambda= m$result[1])
}


mytest=function(snp.now=snpnum[,j]) {
  dat=cbind(dat,snp.now=snp.now)
  SAGE=splines::bs(dat$AGEBL,knots=quantile(dat$AGEBL,c(0.333,0.666)))
  colnames(SAGE)=c("SA1","SA2","SA3","SA4","SA5")
  dat=cbind(dat,SAGE)
  na.sig=rep(1,8)
  #### 1, for all black subjects, transformed
  dat_sub=subset(dat,(dat$dmgrace=="Black"))
  cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
  if (cond) {
    fum_log=hiv~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ts2_aa+I(snp.now*ts2_aa)
    fum_rev=ts2_aa~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
    m1=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
    m2=summary(lm(fum_rev,data=dat_sub))$coef
    
    fum_log=hiv~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ta2_aa+I(snp.now*ta2_aa)
    fum_rev=ta2_aa~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
    m3=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
    m4=summary(lm(fum_rev,data=dat_sub))$coef
    na.sig[1:4]=0
  }
  #### 2, for all Black subjects with exposure >0, transformed
  
  dat_sub=subset(dat,(dat$dmgrace=="Black") & (dat$smoke>0.01))
  cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
  if(cond) {
    fum_log=hiv~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ts2_allc+I(snp.now*ts2_allc)
    fum_rev=ts2_allc~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
    m5=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
    m6=summary(lm(fum_rev,data=dat_sub))$coef
    na.sig[5:6]=0
  }
  
  dat_sub=subset(dat,(dat$dmgrace=="Black") & (dat$alchol>0.01))
  cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
  if(cond) {
    fum_log=hiv~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ta2_allc+I(snp.now*ta2_allc)
    fum_rev=ta2_allc~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
    m7=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
    m8=summary(lm(fum_rev,data=dat_sub))$coef
    na.sig[7:8]=0
  }
  p1=ifelse(na.sig[1],NA,m1[dim(m1)[1],4])
  p2=ifelse(na.sig[2],NA,m2[dim(m2)[1],4])
  p3=ifelse(na.sig[3],NA,m3[dim(m3)[1],4])
  p4=ifelse(na.sig[4],NA,m4[dim(m4)[1],4])
  p5=ifelse(na.sig[5],NA,m5[dim(m5)[1],4])
  p6=ifelse(na.sig[6],NA,m6[dim(m6)[1],4])
  p7=ifelse(na.sig[7],NA,m7[dim(m7)[1],4])
  p8=ifelse(na.sig[8],NA,m8[dim(m8)[1],4])
  
  c(p1,p2,p3,p4,p5,p6,p7,p8)
}






library(foreach)
library(doMC)
library(doRNG)
library("snpStats")
library("car")
library("caret")
library(sandwich)
library(lmtest)
library("GEint")

registerDoMC(16)
registerDoRNG(1)

Sys.time()
setwd("/gpfs/ycga/project/wang_zuoheng/cc2747/VACS/spline200117")
load("/gpfs/ycga/home/cc2747/R_code/covariates.RData")
ethnic=read.csv("/gpfs/ycga/project/wang_zuoheng/ww363/VACs/forCCheng/vacs_selfReport_ethic.csv")[,c(2,5)]
dat=merge(dat,ethnic,by="dnaid",all.x=T)
attach(dat)
dat$alchol= apply(cbind(aBL,aFU1,aFU2,aFU3,aFU4,aFU5),1,mean,na.rm=T)
dat$smoke= apply(cbind(sBL,sFU1,sFU2,sFU3,sFU4,sFU5),1,mean,na.rm=T)
dat$cocaine=apply(cbind(dcq2coke.BL,dcq2coke.FU1,dcq2coke.FU2,dcq2coke.FU3,dcq2coke.FU4,dcq2coke.FU5),1,mean,na.rm=T)
dat$smoke[which(dat$smoke>150)]=NA
dat$dmgrace="Black"

exposure.mat=matrix(NA,ncol=8,nrow=dim(dat)[1])
colnames(exposure.mat)=c("ts1","ts2_all","ts2_aa","ts2_allc","ta1","ta2_all","ta2_aa","ta2_allc")
exposure.mat[,"ts1"]=as.numeric(dat$smoke)
exposure.mat[,"ta1"]=as.numeric(dat$alchol)
exposure.mat[,"ts2_all"]=trans.yj(dat$smoke)
exposure.mat[,"ta2_all"]=trans.yj(dat$alchol)
exposure.mat[which(dat$dmgrace=="Black"),"ts2_aa"]=trans.yj(dat$smoke[which(dat$dmgrace=="Black")])
exposure.mat[which(dat$dmgrace=="Black"),"ta2_aa"]=trans.yj(dat$alchol[which(dat$dmgrace=="Black")])
exposure.mat[which((dat$dmgrace=="Black") & (dat$smoke>0.01)),"ts2_allc"]=trans.yj(dat$smoke[which((dat$dmgrace=="Black") & (dat$smoke>0.01))])
exposure.mat[which((dat$dmgrace=="Black") & (dat$alchol>0.01)),"ta2_allc"]=trans.yj(dat$alchol[which((dat$dmgrace=="Black") & (dat$alchol>0.01))])

dat=cbind(dat,exposure.mat)
i=1
nBS=342

aid <- paste(Sys.getenv("PBS_ARRAYID"),
             Sys.getenv("SLURM_ARRAY_TASK_ID"),
             sep = ""
)
aid = as.numeric(aid)[1]


if (aid < 12) {
  subs= ((aid-1)*29+1):(aid*29)
} else {
  subs=((aid-1)*29+1):342
}




print(paste("aid is",aid))

output <- foreach(i = (subs),.combine = rbind,.errorhandling="remove") %dopar% {
  
  mytest=function(snp.now=snpnum[,j]) {
    dat=cbind(dat,snp.now=snp.now)
    SAGE=splines::bs(dat$AGEBL,knots=quantile(dat$AGEBL,c(0.333,0.666)))
    colnames(SAGE)=c("SA1","SA2","SA3","SA4","SA5")
    dat=cbind(dat,SAGE)
    na.sig=rep(1,8)
    #### 1, for all black subjects, transformed
    dat_sub=subset(dat,(dat$dmgrace=="Black"))
    cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
    if (cond) {
      fum_log=hiv~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ts2_aa+I(snp.now*ts2_aa)
      fum_rev=ts2_aa~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
      m1=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
      m2=summary(lm(fum_rev,data=dat_sub))$coef
      
      fum_log=hiv~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ta2_aa+I(snp.now*ta2_aa)
      fum_rev=ta2_aa~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
      m3=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
      m4=summary(lm(fum_rev,data=dat_sub))$coef
      na.sig[1:4]=0
    }
    #### 2, for all Black subjects with exposure >0, transformed
    
    dat_sub=subset(dat,(dat$dmgrace=="Black") & (dat$smoke>0.01))
    cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
    if(cond) {
      fum_log=hiv~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ts2_allc+I(snp.now*ts2_allc)
      fum_rev=ts2_allc~ta1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
      m5=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
      m6=summary(lm(fum_rev,data=dat_sub))$coef
      na.sig[5:6]=0
    }
    
    dat_sub=subset(dat,(dat$dmgrace=="Black") & (dat$alchol>0.01))
    cond= sum(dat_sub$snp.now,na.rm=T)/(2*length(na.omit(dat_sub$snp.now))) < 0.99
    if(cond) {
      fum_log=hiv~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+ta2_allc+I(snp.now*ta2_allc)
      fum_rev=ta2_allc~ts1+SA1+SA2+SA3+SA4+SA5+dmgsex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+snp.now+hiv+I(snp.now*hiv)
      m7=summary(glm(fum_log,data=dat_sub,family=binomial(link="logit")))$coef
      m8=summary(lm(fum_rev,data=dat_sub))$coef
      na.sig[7:8]=0
    }
    p1=ifelse(na.sig[1],NA,m1[dim(m1)[1],4])
    p2=ifelse(na.sig[2],NA,m2[dim(m2)[1],4])
    p3=ifelse(na.sig[3],NA,m3[dim(m3)[1],4])
    p4=ifelse(na.sig[4],NA,m4[dim(m4)[1],4])
    p5=ifelse(na.sig[5],NA,m5[dim(m5)[1],4])
    p6=ifelse(na.sig[6],NA,m6[dim(m6)[1],4])
    p7=ifelse(na.sig[7],NA,m7[dim(m7)[1],4])
    p8=ifelse(na.sig[8],NA,m8[dim(m8)[1],4])
    
    c(p1,p2,p3,p4,p5,p6,p7,p8)
  }
  
  if (i<342) {
    loc=((i-1)*50000+1):((i-1)*50000+50000)
  } else {
    loc=((i-1)*50000+1):17092657
  }  
  snp=snpStats::read.plink(bed="/ycga-gpfs/project/fas/wang_zuoheng/ww363/VACs/genoQC/hwe6_noduplicate_nooutlier_QC_2.bed", 
                           bim="/ycga-gpfs/project/fas/wang_zuoheng/ww363/VACs/genoQC/hwe6_noduplicate_nooutlier_QC_2.bim", 
                           fam="/ycga-gpfs/project/fas/wang_zuoheng/ww363/VACs/genoQC/hwe6_noduplicate_nooutlier_QC_2.fam",
                           select.snps=loc)
  # The genotype matrix
  snpmatrix <- snp$genotypes
  # creat SNP summary statistics
  colsum = snpStats::col.summary(snpmatrix)
  ## setting thresholds
  ## MAF > minor
  minor=0.01
  # call rate > call
  call=0.95
  ## filter on MAF and call rate
  use=with(colsum,(!is.na(MAF) & MAF>minor) & Call.rate >=call)
  ## remove NAs
  use[is.na(use)]=FALSE
  
  ## subset genotype and SNP summary data
  snpmatrix=snpmatrix[,which(use==1)]
  manhattan=snp$map[which(use==1),c(1,2,4)]
  colsum=colsum[which(use==1),]
  #print(snpmatrix)
  snpnum=as(snpmatrix,"numeric")
  rm(snpmatrix)
  rownames(snpnum)=snp$fam$member
  rownames(dat)=dat$dnaid
  common.id=intersect(as.integer(snp$fam$member),as.integer(dat$dnaid))
  snpnum=snpnum[as.integer(rownames(snpnum)) %in% common.id,]
  dat=dat[as.integer(rownames(dat)) %in% common.id,]
  snpnum=snpnum[order(as.integer(rownames(snpnum))),]
  
  rm(snp)
  ## linear model
  n=dim(snpnum)[2]
  pvalue.mat=matrix(0.5,ncol=8,nrow=n)
  for (j in (1:n)) {
    pvalue.mat[j,]=mytest(snp.now = snpnum[,j])
    
  }
  manhattan=cbind(manhattan,pvalue.mat)
  save(manhattan, file = paste("The ",i,"-th manhattan data trans.RData",sep=""))
  rm(manhattan)
  rm(snpnum)
  print(i)
  #gc()
  return(1)
}
Sys.time()
















