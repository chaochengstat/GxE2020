setwd("/gpfs/ycga/project/wang_zuoheng/cc2747/VACS/simulation/P_05P/")
library(foreach)
library(doMC)
library(doRNG)
library(sn)
test.fun = function(rho=1,P_1000=0.05, alpha_g=log(1.5), alpha_x=log(1.5), alpha_gx=log(1.2),
                    Cor_XZ1=0.1, Cor_XZ2=0.1, Cor_XG=0.3) {
  #rho=1;P_1000=0.05; alpha_g=log(1.5); alpha_x=log(1.5); alpha_gx=log(1.2)
  #Cor_XZ1=0.2; Cor_XZ2=0.1; Cor_XG=0.3
  ## NO of individuals in the whole cohort study
  N  = 30000
  ## Generate Z
  #Z1=rnorm(N,sd=2)
  #plot(Z1,1/(1+exp(Z1)))
  
  Z1=rnorm(N,sd=1)
  Z2= as.numeric(runif(N) < 0.3)
  ## Generate G
  a1 = log(0.3/0.7)+log(1.2)*Z1+log(1.2)*Z2
  a2 = log(0.1/0.9)+log(1.5)*Z1+log(1.5)*Z2
  P_G1=exp(a1)/(1+exp(a1)+exp(a2))
  P_G2=exp(a2)/(1+exp(a1)+exp(a2))
  P_G0=1-P_G1-P_G2
  Gr = runif(N)
  G = (Gr<P_G2)*1+(Gr<(P_G2+P_G1))*1
  ## Generate X
  alpha_0=log(P_1000/(1-P_1000))
  alpha_1=log(1.2)
  alpha_2=log(1.2)
  alpha_g=alpha_g
  alpha_x=alpha_x
  alpha_gx=alpha_gx
  
  Cor_XZ1= Cor_XZ1
  Cor_XZ2= Cor_XZ2
  Cor_XG = Cor_XG
  beta_0 = 0
  beta_1 = Cor_XZ1/sd(Z1+Z1^2)
  beta_2 = Cor_XZ2/sd(Z2)
  beta_g = Cor_XG/sd(G)
  beta_d = alpha_x
  beta_gd = alpha_gx
  SZ1=Z1+Z1^2
  d1=0.5*(alpha_x+alpha_gx*G)^2
  d2=alpha_0 + alpha_1*SZ1 + alpha_2*Z2 + alpha_g*G
  d3=alpha_x + alpha_gx*G
  d4=beta_0 + beta_1*SZ1 + beta_2*Z2 + beta_g*G
  P_D = 1/(1+exp(-(d1+d2+d3*d4)))
  Ds = as.numeric(runif(N) < P_D)
  Mean_X0 = d4
  Mean_X1 = d4+beta_d+beta_gd*G
  Mean_X  = Mean_X0*(1-Ds) + Mean_X1 * Ds
  X = Mean_X+rnorm(N)
  if (rho==1) {
    Xm=X
  } else {
    Xm=X+rnorm(length(X),mean=0,sd=sqrt(1/rho - 1))
  }
  
  ## Generate D
  P_D = 1/(1+exp(-(d2+d3*X)))
  D = as.numeric(runif(N) < P_D)
  
  # cohort study
  data.cohort=data.frame(Z1=Z1,Z2=Z2,X=X,G=G,D=D,Xm=Xm)
  # case control study
  n=1000
  if(sum(D)<n) {
    n=sum(D)
  }
  case.index=sample(which(D==1),n,replace=F)
  control.index=sample(which(D==0),n,replace=F)
  #data.casecontrol = data.cohort[c(case.index,control.index),]
  
  # log test
  index=c(case.index,control.index)
  Z1=Z1[index]
  Z2=Z2[index]
  G =G[index]
  D =D[index]
  Xm=Xm[index]
  SZ1=splines::bs(Z1,knots=quantile(Z1,c(0.333,0.666)))
  #SZ2=splines::bs(Z2,knots=quantile(Z2,0.333,0.666))
  S1 = glm(D~ SZ1+Z2 + Xm + G + I(Xm*G) , family=binomial(link="logit"))
  r=dim(summary(S1)$coefficients)[1]
  log.test=summary(S1)$coefficients[r,c(3,4)];log.test[1]=log.test[1]^2
  # rev test
  S2 = lm(Xm ~ SZ1 + Z2 + D + G + I(D*G))
  r=dim(summary(S2)$coefficients)[1]
  rev.test=summary(S2)$coefficients[r,c(3,4)];rev.test[1]=rev.test[1]^2
  output=data.frame(log.test=log.test,rev.test=rev.test)
  if ( log.test[1]>1000 ) {
    output=data.frame(log.test=c(1,0.05),rev.test=c(1,0.05))
  }
  rownames(output)=c("Chi-Squared","P-Value")
  output
}


test.fun.main=function(rep=100,rho=0.5,P_1000=0.05, alpha_g=log(1.5), alpha_x=log(1.5), alpha_gx=log(1.0),
                       Cor_XZ1=0.1, Cor_XZ2=0.0, Cor_XG=0.3) {
  #P_log=rep(NA,rep)
  #P_rev=rep(NA,rep)
  #Chi_log=Chi_rev=rep(NA,rep)
  #rep=100;rho=0.5;P_1000=0.05; alpha_g=log(1.5); alpha_x=log(1.5); alpha_gx=log(1.0);
  #Cor_XZ1=0.1; Cor_XZ2=0.0; Cor_XG=0.3
  cl <- 15
  registerDoMC(cl)
  registerDoRNG(1)
  res=foreach(i = (1:rep),.combine = rbind,.errorhandling="remove") %dopar% {
    
    test.fun = function(rho=1,P_1000=0.05, alpha_g=log(1.5), alpha_x=log(1.5), alpha_gx=log(1.2),
                        Cor_XZ1=0.1, Cor_XZ2=0.1, Cor_XG=0.3) {
      #rho=1;P_1000=0.05; alpha_g=log(1.5); alpha_x=log(1.5); alpha_gx=log(1.2)
      #Cor_XZ1=0.2; Cor_XZ2=0.1; Cor_XG=0.3
      ## NO of individuals in the whole cohort study
      N  = 30000
      ## Generate Z
      #Z1=rnorm(N,sd=2)
      #plot(Z1,1/(1+exp(Z1)))
      
      Z1=rnorm(N,sd=1)
      Z2= as.numeric(runif(N) < 0.3)
      ## Generate G
      a1 = log(0.3/0.7)+log(1.2)*Z1+log(1.2)*Z2
      a2 = log(0.1/0.9)+log(1.5)*Z1+log(1.5)*Z2
      P_G1=exp(a1)/(1+exp(a1)+exp(a2))
      P_G2=exp(a2)/(1+exp(a1)+exp(a2))
      P_G0=1-P_G1-P_G2
      Gr = runif(N)
      G = (Gr<P_G2)*1+(Gr<(P_G2+P_G1))*1
      ## Generate X
      alpha_0=log(P_1000/(1-P_1000))
      alpha_1=log(1.2)
      alpha_2=log(1.2)
      alpha_g=alpha_g
      alpha_x=alpha_x
      alpha_gx=alpha_gx
      
      Cor_XZ1= Cor_XZ1
      Cor_XZ2= Cor_XZ2
      Cor_XG = Cor_XG
      beta_0 = 0
      beta_1 = Cor_XZ1/sd(Z1+Z1^2)
      beta_2 = Cor_XZ2/sd(Z2)
      beta_g = Cor_XG/sd(G)
      beta_d = alpha_x
      beta_gd = alpha_gx
      SZ1=Z1+Z1^2
      d1=0.5*(alpha_x+alpha_gx*G)^2
      d2=alpha_0 + alpha_1*SZ1 + alpha_2*Z2 + alpha_g*G
      d3=alpha_x + alpha_gx*G
      d4=beta_0 + beta_1*SZ1 + beta_2*Z2 + beta_g*G
      P_D = 1/(1+exp(-(d1+d2+d3*d4)))
      Ds = as.numeric(runif(N) < P_D)
  Mean_X0 = d4
  Mean_X1 = d4+beta_d+beta_gd*G
  Mean_X  = Mean_X0*(1-Ds) + Mean_X1 * Ds
  X = rnorm(N,mean=Mean_X,sd=1)
      if (rho==1) {
        Xm=X
      } else {
        Xm=X+rnorm(length(X),mean=0,sd=sqrt(1/rho - 1))
      }
      
      ## Generate D
      P_D = 1/(1+exp(-(d2+d3*X)))
      D = as.numeric(runif(N) < P_D)
      
      # cohort study
      data.cohort=data.frame(Z1=Z1,Z2=Z2,X=X,G=G,D=D,Xm=Xm)
      # case control study
      n=1000
      if(sum(D)<n) {
        n=sum(D)
      }
      case.index=sample(which(D==1),n,replace=F)
      control.index=sample(which(D==0),n,replace=F)
      #data.casecontrol = data.cohort[c(case.index,control.index),]
      
      # log test
      index=c(case.index,control.index)
      Z1=Z1[index]
      Z2=Z2[index]
      G =G[index]
      D =D[index]
      Xm=Xm[index]
      SZ1=splines::bs(Z1,knots=quantile(Z1,c(0.333,0.666)))
      #SZ2=splines::bs(Z2,knots=quantile(Z2,0.333,0.666))
      S1 = glm(D~ SZ1+Z2 + Xm + G + I(Xm*G) , family=binomial(link="logit"))
      r=dim(summary(S1)$coefficients)[1]
      log.test=summary(S1)$coefficients[r,c(3,4)];log.test[1]=log.test[1]^2
      # rev test
      S2 = lm(Xm ~ SZ1 + Z2 + D + G + I(D*G))
      r=dim(summary(S2)$coefficients)[1]
      rev.test=summary(S2)$coefficients[r,c(3,4)];rev.test[1]=rev.test[1]^2
      output=data.frame(log.test=log.test,rev.test=rev.test)
      if ( log.test[1]>1000 ) {
        output=data.frame(log.test=c(1,0.05),rev.test=c(1,0.05))
      }
      rownames(output)=c("Chi-Squared","P-Value")
      output
    }
    
    
    res=test.fun(rho=rho,P_1000=P_1000, alpha_g=alpha_g, alpha_x=alpha_x, alpha_gx=alpha_gx,
                 Cor_XZ1=Cor_XZ1, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
    P_log=res[2,1];P_rev=res[2,2]
    Chi_log=res[1,1];Chi_rev=res[1,2]
    c(P_log,P_rev,Chi_log,Chi_rev)
  }
  return(c(IE_log=mean(res[,1] < 0.05),IE_rev=mean(res[,2] < 0.05),Chi2_ratio=mean(res[,4])/mean(res[,3])))
}


all.comb=function(x) {
  #x=c(2,2,2,3)
  output=matrix(0,ncol=length(x),nrow=cumprod(x)[length(x)])
  z=x
  for (i in (1:length(x))) {
    z=z[-1]
    output[,i]=rep(1:x[i],each=cumprod(z)[length(z)])
  }
  output
}

print(Sys.time())
options(warn=-1)
rep=500000#000
P_1000_list  = c(0.05,0.2,0.5)
alpha_g_list = c(log(1.1),log(1.5))
#alpha_x_list = c(log(1.1),log(1.5))
Cor_XZ2_list = c(0.01,0.2)
Cor_XG_list  = c(0.01,0.2)




###########################################################
# Index of Task
###########################################################
aid <- paste(Sys.getenv("PBS_ARRAYID"),
             Sys.getenv("SLURM_ARRAY_TASK_ID"),
             sep = ""
)
aid = as.numeric(aid)[1]

if (aid < 8) {
  subs=((aid-1)*3+1):(aid*3)
} else {
  subs=22:24
}



#beta_X_list=c(log(1),log(1.5))
allcomb=all.comb(c(3,2,2,2))
res=matrix(NA,ncol=8,nrow=dim(allcomb)[1])
par.mat=matrix(NA,ncol=dim(allcomb)[2],nrow=dim(allcomb)[1])
for (i in (subs)) {
  parm.index=allcomb[i,]
  P_1000=P_1000_list[parm.index[1]]
  alpha_g=alpha_g_list[parm.index[2]]
  Cor_XZ2 = Cor_XZ2_list[parm.index[3]]
  Cor_XG=Cor_XG_list[parm.index[4]]
  test1=test.fun.main(rep=rep,rho=1,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.1), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test2=test.fun.main(rep=rep,rho=0.75,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.1), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test3=test.fun.main(rep=rep,rho=0.5,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.1), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test4=test.fun.main(rep=rep,rho=0.25,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.1), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  
  test5=test.fun.main(rep=rep,rho=1,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.5), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test6=test.fun.main(rep=rep,rho=0.75,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.5), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test7=test.fun.main(rep=rep,rho=0.5,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.5), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  test8=test.fun.main(rep=rep,rho=0.25,P_1000=P_1000, alpha_g=alpha_g, alpha_x=log(1.5), alpha_gx=log(1.5),
                      Cor_XZ1=0.01, Cor_XZ2=Cor_XZ2, Cor_XG=Cor_XG)
  res[i,]=c(test1[3],test2[3],test3[3],test4[3],test5[3],test6[3],test7[3],test8[3])
  par.mat[i,]=c(P_1000,alpha_g,Cor_XZ2,Cor_XG)
  print(paste("NO.",i,"and",Sys.time()))
  print(res[i,])
}
file=paste("Power_skew_",aid,".csv",sep="")
write.csv(res,file)




